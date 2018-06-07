
from distutils.core import setup

setup(name='ampel-contrib-hu',
      version='0.1',
      packages=['ampel.contrib.hu',
                'ampel.contrib.hu.examples.t0',
                'ampel.contrib.hu.examples.t2',
                'ampel.contrib.hu.t0',
                'ampel.contrib.hu.t2',
                'ampel.contrib.hu.t3'],
      package_data = {'': ['*.json']},
      entry_points = {
          'ampel.channels' : [
              'nuclear_transient_charlotte = ampel.contrib.ztfbh.config:charlotte',
              'nuclear_transient_sjoert = ampel.contrib.ztfbh.config:sjoert',
          ],
          'ampel.pipeline.t0' : [
              'DecentFilter = ampel.contrib.hu.t0.DecentFilter:DecentFilter',
              'LensedTransientFilter = ampel.contrib.hu.t0.LensedTransientFilter:LensedTransientFilter',
              'NeutrinoFilter = ampel.contrib.hu.t0.NeutrinoFilter:NeutrinoFilter',
              'NoFilter = ampel.contrib.hu.t0.NoFilter:NoFilter',
              'RandFilter = ampel.contrib.hu.t0.RandFilter:RandFilter',
              'SEDmTargetFilter = ampel.contrib.hu.t0.SEDmTargetFilter:SEDmTargetFilter',
              'SNFilter = ampel.contrib.hu.t0.SNFilter:SNFilter',
          ]
      }
)
