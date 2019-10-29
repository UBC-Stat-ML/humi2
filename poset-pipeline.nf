#!/usr/bin/env nextflow

import java.util.stream.Collectors



pwd=new File(".").getAbsolutePath()

datasets = Channel.fromPath( 'data/*', type: 'dir' ).filter{ 
  result = !it.getName().startsWith(".") 
  if (result) println("Queuing data " + it.toString())
  return result
}.map{it.getName().toString()}

params.model = "BNB LocalLambdaMixBNB MixBNB MixNB MixYS NB Poi YS"

params.nScans = 1000
params.nInitParticles = 10 // increase this if model initialization fails (can happen in complex mixture models with vague priors)
params.nTargets = "INF" // use this to do inference on a subset of targets (e.g. for dry runs)
params.nExperiments = "INF"

deliverableDir = 'deliverables/' + workflow.scriptName.replace('.nf','') + "_" + params.nScans + "_" + params.nInitParticles + "_" + params.nTargets + "_" + params.nExperiments + "/"
runsDir = deliverableDir + "runs" 
posetsDir = deliverableDir + "posets" 

models = Arrays.asList(params.model.split("\\s+")).stream().map{
  result = "humi.models." + it
  println("Queuing model " + result)
  return result
}.collect(Collectors.toList())

process buildCode {
  cache true 
  input:
    val gitRepoName from 'nowellpack'
    val gitUser from 'UBC-Stat-ML'
    val codeRevision from '5b84c0aa2255c4cb932517bf16b52656e0a1eadc'
    val snapshotPath from "${System.getProperty('user.home')}/w/nowellpack"
  output:
    file 'code' into code
  script:
    template 'buildRepo.sh' 
}

process run {

  // uncomment to use on cluster:
  // cpus 1
  // executor 'sge'
  // memory '5 GB'
  // time '50h'

  input:
    file code
    each model from models
    each dataset from datasets
  output:
    file "${dataset}_${model}" into runs
  publishDir runsDir, mode: 'link'
  """
  ln -s $pwd/data/$dataset/final.csv ${dataset}.csv  # o.w. spark crash on it
  java -cp code/lib/\\* -Xmx1g $model   \
           --model.initialPopCounts.dataSource $pwd/data/$dataset/initial.csv \
           --model.initialPopCounts.name counts  \
           --model.data.source ${dataset}.csv \
           --model.data.genes.name gene     \
           --model.data.targets.name sgRNA     \
           --model.data.targets.maxSize $params.nTargets \
           --model.data.experiments.name dataset     \
           --model.data.experiments.maxSize ${params.nExperiments} \
           --model.data.histograms.name histogram     \
           --engine.nScans $params.nScans   \
           --engine.nChains 1 \
           --engine.nPassesPerScan 1     \
           --engine.nThreads Fixed     \
           --engine.nThreads.number 1 \
           --engine.scmInit.nParticles $params.nInitParticles \
           --engine.scmInit.temperatureSchedule.threshold 0.6 \
           --engine.scmInit.nThreads Fixed \
           --engine.scmInit.nThreads.number 1 \
           --postProcessor humi.HumiPostProcessor \
           --postProcessor.data.targets.name sgRNA \
           --postProcessor.data.genes.name gene \
           --postProcessor.data.experiments.name dataset \
           --postProcessor.data.histograms.name histogram \
           --postProcessor.runPxviz false
  mv results/all/`ls results/all` ${dataset}_${model}
  """
}

process poset {
  input:
    file code
    file run from runs
  output:
    file "*.dot"
  publishDir posetsDir, mode: 'copy', overwrite: true
  """
  java -cp code/lib/\\* -Xmx1g humi.posets.Intervals2Poset \
    --intervalsCSVFile $run/estimates.csv
  mv results/latest/hasse.dot ${run}.dot
  """
}

process analysisCode {
  input:
    val gitRepoName from 'nedry'
    val gitUser from 'alexandrebouchard'
    val codeRevision from 'cf1a17574f19f22c4caf6878669df921df27c868'
    val snapshotPath from "${System.getProperty('user.home')}/w/nedry"
  output:
    file 'code' into analysisCode
  script:
    template 'buildRepo.sh'
}



process summarizePipeline {
  cache false 
  output:
      file 'pipeline-info.txt'
  publishDir deliverableDir, mode: 'copy', overwrite: true
  """
  echo 'scriptName: $workflow.scriptName' >> pipeline-info.txt
  echo 'start: $workflow.start' >> pipeline-info.txt
  echo 'runName: $workflow.runName' >> pipeline-info.txt
  echo 'nextflow.version: $workflow.nextflow.version' >> pipeline-info.txt
  """
}
