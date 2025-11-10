mamba activate dualSeq

python dualSeq/dualSeq.py --data sample_overview.csv --pipeline-config pipeline.yaml --outdir results/ --profile snakemake_profile/ -t 8 -j 16 --use-conda --latency-wait 60

