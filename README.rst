Contributed Ampel units from HU/DESY group
==========================================

Demo install instructions:
==========================

Create environment with python 3.10+ / poetry. Then run:
- `git clone https://github.com/AmpelProject/Ampel-HU-astro.git`
- `cd Ampel-HU-astro/`
- `poetry install -E "ztf sncosmo extcats notebook"`
- `cd notebooks`
- `poetry run jupyter notebook`

This will allow a number of _Demo_ notebooks to be run. Note that most of them
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
