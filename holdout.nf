#!/usr/bin/env nextflow

import java.util.stream.Collectors



pwd=new File(".").getAbsolutePath()

params.dataset = "X0569_422234x"
params.maxNExperiments = 4
params.lastDataset = "X0569_4222344_clones"

expSizes = 1..params.maxNExperiments

params.model = "BNB LocalLambdaMixBNB MixBNB MixNB MixYS NB Poi YS"

params.nScans = 1000
params.nInitParticles = 10 // increase this if model initialization fails (can happen in complex mixture models with vague priors)
params.nTargets = "INF" // use this to do inference on a subset of targets (e.g. for dry runs)

deliverableDir = 'deliverables/' + workflow.scriptName.replace('.nf','') + "_" + params.nScans + "_" + params.nInitParticles + "_" + params.nTargets + "/" + params.dataset  + "/"
runsDir = deliverableDir + "runs" 

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
    each expSize from expSizes
  output:
    file "${model}_${expSize}" into runs
  publishDir runsDir, mode: 'link'
"""
  java -cp code/lib/\\* -Xmx5g $model   \
           --model.initialPopCounts.dataSource $pwd/data/$params.dataset/initial.csv \
           --model.initialPopCounts.name counts  \
           --model.data.source $pwd/data/$params.dataset/final.csv \
           --model.data.genes.name gene     \
           --model.data.targets.name sgRNA     \
           --model.data.targets.maxSize $params.nTargets \
           --model.data.experiments.name dataset     \
           --model.data.experiments.maxSize $expSize \
           --model.data.histograms.name histogram     \
           --engine.nScans $params.nScans   \
           --engine.nChains 1 \
           --engine.nPassesPerScan 1     \
           --engine.nThreads Fixed     \
           --engine.nThreads.number 1 \
           --engine.scmInit.nParticles $params.nInitParticles \
           --engine.scmInit.temperatureSchedule.threshold 0.9 \
           --engine.scmInit.nThreads Fixed \
           --engine.scmInit.nThreads.number 1 \
           --postProcessor humi.HumiPostProcessor \
           --postProcessor.data.targets.name sgRNA \
           --postProcessor.data.genes.name gene \
           --postProcessor.data.experiments.name dataset \
           --postProcessor.data.histograms.name histogram \
           --postProcessor.runPxviz false
  mv results/all/`ls results/all` ${model}_${expSize}
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


process aggregate {
  input:
    file analysisCode
    file 'exec_*' from runs.toList()
  output:
    file 'results/latest/aggregated' into aggregated
  """
  code/bin/aggregate \
    --dataPathInEachExecFolder gof.csv \
    --keys \
      model.data.experiments.maxSize as nExperiments \
      model \
        from arguments.tsv
  """
}

process plot {
  input:
    file aggregated
    env SPARK_HOME from "${System.getProperty('user.home')}/bin/spark-2.1.0-bin-hadoop2.7"
  output:
    file '*.csv'
    file '*.pdf'
  publishDir deliverableDir, mode: 'copy', overwrite: true
  afterScript 'rm -r metastore_db; rm derby.log'
  """
  #!/usr/bin/env Rscript
  require("ggplot2")
  require("stringr")
  library(SparkR, lib.loc = c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib")))
  sparkR.session(master = "local[*]", sparkConfig = list(spark.driver.memory = "4g"))

  data <- read.df("$aggregated", "csv", header="true", inferSchema="true")
  data <- collect(data)
  
  require("dplyr")
  
  data\$model <- str_replace_all(data\$model, "[\$].*", "")
  data\$model <- str_replace_all(data\$model, "humi[.]models[.]", "")
  
  write.csv(data, file="gd-data.csv")
  
  data <- data %>%
    filter(gofStatistic != "visibleCloneNumbers") %>%
    filter(referenceDataset == "${params.lastDataset}")
    
  p <- ggplot(data, aes(x=nExperiments, y=actualCoverage, colour = model, group = model)) + 
    geom_line() +
    facet_grid(gofStatistic ~ ., scales="free") + 
    geom_hline(yintercept=0.95) + 
    theme_bw() 

  ggsave(plot=p, file="generalization-coverage.pdf")
 

  p <- ggplot(data, aes(x=nExperiments, y=width, , colour = model, group = model)) + 
    geom_line() +
    facet_grid(gofStatistic ~ ., scales="free") + 
    theme_bw() 

  ggsave(plot=p, file="generalization-width.pdf")  
  """
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
