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
  executor 'local'
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

  time '1h'

  input:
    file code
    each model from models
    each dataset from datasets
    each vague from 0.1, 0.01, 0.001
  output:
    file "results/latest" into runs_
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
           --model.vagueRate $vague \
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
  
  """
}

runs_.into {
  runs
  to_aggregate
}

process analysisCode {
  executor 'local'
  cache true 
  input:
    val gitRepoName from 'nedry'
    val gitUser from 'UBC-Stat-ML'
    val codeRevision from 'a9abcc40abcfb285588cc4c312d8ecc0bbdad06e'
    val snapshotPath from "${System.getProperty('user.home')}/w/nedry"
  output:
    file 'code' into analysisCode
  script:
    template 'buildRepo.sh'
}

process aggregate {
  echo true
  scratch false
  input:
    file analysisCode
    file 'exec_*' from to_aggregate.toList()
  output:
    file 'results/latest/estimates.csv' into aggregated
  """
  java -Xmx5g -cp code/lib/\\*  flows.Aggregate \
    --experimentConfigs.resultsHTMLPage false \
    --dataPathInEachExecFolder estimates.csv \
    --keys \
      model.data.source as dataset \
      model.vagueRate as priorHyperParameter \
           from arguments.tsv
  """
}

process plot {
  scratch false
  container 'cgrlab/tidyverse'
  input:
    file aggregated
  output:
    file '*.pdf'
  publishDir deliverableDir, mode: 'copy', overwrite: true
  """
  #!/usr/bin/env Rscript
  require("ggplot2")
  require("dplyr")
  
  data <- read.csv("$aggregated")
  
  
  cols = rainbow(200, s=.6, v=.9)[sample(1:200,200)]
  p <- ggplot(data, aes(x = factor(priorHyperParameter), y = logRatio, colour = factor(sgrna))) + 
    coord_flip() + 
    geom_errorbar(aes(ymin=logRatioLeftBound, ymax=logRatioRightBound)) +
    geom_point() + 
    facet_grid(factor(gene) + factor(sgrna)  ~ dataset) +
    theme_bw() + 
    xlab("Gene") + 
    ylab("log(ratio)") + 
    ggtitle("Ratio of clone sizes relative to controls", subtitle = "Bayesian hierarchical model credible intervals") + 
    scale_colour_manual(values=cols) + 
    geom_hline(yintercept=0) + 
    theme(legend.position="none") 
  ggsave(plot = p, filename = "intervals-sensitivity.pdf", height = 100, limitsize = FALSE)

  """
}

