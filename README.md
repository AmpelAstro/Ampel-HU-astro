<img align="left" src="https://user-images.githubusercontent.com/17532220/213287034-0209aa19-f8a1-418f-a325-7472510542cb.png" width="150" height="150"/>
<br>

# AMPEL-HU-astro
<br><br>


Contributed Ampel units from HU/DESY group
==========================================

Demo install instructions:
==========================

1. [Install poetry](https://python-poetry.org/docs/#installation). If you install poetry with conda, be sure to install it in its own environment, e.g. `conda create -n poetry`.
2. `git clone https://github.com/AmpelProject/Ampel-HU-astro.git; cd Ampel-HU-astro`
3. Check your virtualenv setup with `poetry env info` (or `conda run -n poetry poetry env info` if using conda). The output should include:
   ```
   Virtualenv
   Python:         3.10.x
   ```
   If not, point poetry at an installation of Python 3.10 with (`conda run -n poetry)` `poetry env use PATH_TO_PYTHON_310`
4. (`conda run -n poetry`) `poetry install -E "ztf sncosmo extcats notebook"`
5. `cd notebooks`
6. (`conda run -n poetry`) `poetry run jupyter notebook`

This will allow a number of Demo / access / development notebooks to be run. Note that most of them
requires an access token if data is to be retrieved.

Contains as of Nov 2022:
========================

T0
--
* SimpleDecentFilter
* LensedTransientFilter
* NoFilter
* RandFilter
* SEDmTargetFilter
* SimpleDecentFilter
* ToOFilter
* TransientInClusterFilter
* TransientInEllipticalFilter
* XShooterFilter
* RcfFilter
* RedshiftCatalogFilter

T2
--
* T2PanStarrThumbPrint
* T2PhaseLimit
* T2PS1ThumbExtCat
* T2PS1ThumbNedSNCosmo
* T2PS1ThumbNedTap
* T2LCQuality
* T2BrightSNProb
* T2TNSEval
* T2InfantCatalogEval
* T2RunSncosmo
* T2CatalogMatchLocal
* T2DigestRedshifts
* T2RunPossis
* T2RunTDE
* T2RunParsnip
* T2RunSnoopy
* T2MatchBTS
* T2NedTap
* T2NedSNCosmo
* T2PropagateStockInfo
* T2GetLensSNParameters
* T2LSPhotoZTap
* T2ElasticcRedshiftSampler
* T2TabulatorRiseDecline
* T2XgbClassifier
* T2ElasticcReport
* T2FastDecliner

T3
--
* TransientInfoPrinter
* TransientViewDumper
* ChannelSummaryPublisher
* SlackSummaryPublisher
* RapidBase
* RapidSedm
* RapidLco
* TNSTalker
* TNSMirrorUpdater
* TransientTablePublisher
* HealpixCorrPlotter
* PlotLightcurveSample
* ElasticcClassPublisher
* VOEventPublisher
