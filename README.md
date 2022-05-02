# Branch containing specific tutorial/demo notebooks based on 0.8.2

## Install guidelines

On Mac (M1), additional preparation is needed:
- We need to install Rust and a Kafka dependency with ```brew install rustup librdkafka```, then execute
- ```rustup-init```
- After that, add the following lines to your ~/.zshrc:
    - ```export C_INCLUDE_PATH=/opt/homebrew/Cellar/librdkafka/1.8.2/include``` and
    - ```export LIBRARY_PATH=/opt/homebrew/Cellar/librdkafka/1.8.2/lib``` 
- Now read the modified .zshrc by ```source ~/.zshrc```, then
- execute ```pip3 install confluent_kafka``` and you are good to go.

Sample instructions for creating a conda environment

- ```conda create -n ampelTutorial python=3.9```
- ```conda activate ampelTutorial```
- ```pip3 install ipython jupyter extcats sfdmap iminuit sncosmo light-curve-python```
- ```pip3 install git+https://github.com/AmpelProject/Ampel-ipython.git```
- ```pip3 install git+https://github.com/AmpelProject/Ampel-ZTF.git@dev/v0.8.2 ```
- ```git clone https://github.com/AmpelProject/Ampel-HU-astro.git```
- ```cd Ampel-HU-astro/```
- ```git checkout AmpelTutorial```
- ```pip install -e .```
- ```cd ..```
- ```ampel config build -out ampel_conf.yaml >& conf.log```

The last command will create an Ampel yaml configuration file, which will be required when running any AMPEL context units.

## Access tokens
Access data either from the ZTF alert archive or from the live AMPEL instance requires an access token. This is created through a github auth - contact maintainers to get added to the project.
