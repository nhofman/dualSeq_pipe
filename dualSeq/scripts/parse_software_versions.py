import yaml
import sys
from pathlib import Path

envs_dir = sys.argv[1]
software_list = []
yaml_list = Path(envs_dir).glob("*.yaml")
file = open('multiqc_config.yaml', 'a')
file.write('software_versions:\n')
for conda_env in yaml_list:
    conda = yaml.safe_load(conda_env.open('r'))
    for entry in conda['dependencies']:
        if "=" in entry:
            entry_split = entry.split("=")
            tmp = str(entry_split[0]) + ': "' + str(entry_split[1] + '"')
            file.write(f'  {tmp}\n')
file.close()
