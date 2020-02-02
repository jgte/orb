# YAMLMatlab

This imports YAML files into MATLAB. I found this on the Google Code Archive. I've copied it here so that it can be easily be submoduled into other projects.

## Usage

1. Add the path to this cloned repository: `addpath('path/to/this/repo')`
1. Reading: `YamlStruct = ReadYaml('filename.yaml')`
1. Writing:
    * Create struct: `x.name = 'Martin'`
    * Create YAML: `WriteYaml('filename.yaml',x)`