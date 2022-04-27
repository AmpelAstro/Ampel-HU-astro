Branch containing specific tutorial/demo notebooks based on 0.8.2

Install guidelines
==================

Sample instructions for creating a conda environment

- pip install ipython 
- pip install jupyter 
- pip install extcats
- pip install sfdmap
- pip install iminuit
- pip install sncosmo
- pip install light-curve-python
- conda create -n ampelTutorial python=3.9
- conda activate ampelTutorial
- git clone https://github.com/AmpelProject/Ampel-HU-astro.git 
- cd Ampel-HU-astro/
- git checkout AmpelTutorial
- cd ..
- git clone https://github.com/AmpelProject/Ampel-ipython.git
- cd Ampel-ipython/
- pip install -e .
- ampel config build -out ampel_conf.yaml >& ampel_conf.log 

The last command will create an Ampel yaml configuration file, which will be required when running any AMPEL context units.

Access tokens
=============

Access data either from the ZTF alert archive or from the live AMPEL instance requires an access token. This is created through a github auth - contact maintainers to get added to the project.
