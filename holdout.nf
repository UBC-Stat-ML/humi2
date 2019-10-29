#!/usr/bin/env nextflow

deliverableDir = 'deliverables/' + workflow.scriptName.replace('.nf','')

pwd=new File(".").getAbsolutePath()

params.finalCounts
params.initialCounts
params.maxNExperiments
params.lastDataset

expSizes = 1..params.maxNExperiments

finalCounts = new File(params.finalCounts)
initialCounts = new File(params.initialCounts)

process buildCode {
  cache true 
  input:
    val gitRepoName from 'nowellpack'
    val gitUser from 'UBC-Stat-ML'
    val codeRevision from '1ab29b6288358aa6c0004f830313b3aab2668d5f'
    val snapshotPath from "${System.getProperty('user.home')}/w/nowellpack"
  output:
    file 'code' into code
  script:
    template 'buildRepo.sh' 
}

process run {
  input:
    file code
    each expSize from expSizes
    each model from 'humi.models.NB', 'humi.models.MixBNB', 'humi.models.BNB', 'humi.models.Poi', 'humi.models.MixNB', 'humi.models.GlobalLambdaMixBNB'
  output:
    file 'results/latest' into runs
"""
  java -cp code/lib/\\* -Xmx5g $model   \
           --model.initialPopCounts.dataSource $initialCounts \
           --model.initialPopCounts.name counts  \
           --model.data.source $finalCounts \
           --model.data.genes.name gene     \
           --model.data.targets.name sgRNA     \
           --model.data.targets.maxSize INF \
           --model.data.experiments.name dataset     \
           --model.data.experiments.maxSize $expSize \
           --model.data.histograms.name histogram     \
           --engine.nScans 1000   \
           --engine.nChains 1 \
           --engine.nPassesPerScan 1     \
           --engine.nThreads Fixed     \
           --engine.nThreads.number 1 \
           --engine.scmInit.nParticles 10 \
           --engine.scmInit.temperatureSchedule.threshold 0.9 \
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
