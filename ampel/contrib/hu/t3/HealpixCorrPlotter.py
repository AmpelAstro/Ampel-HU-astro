#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:               Ampel-contrib-HU/ampel/contrib/hu/t3/HealpixCorrPlotter.py
# License:            BSD-3-Clause
# Author:             jn <jnordin@physik.hu-berlin.de>
# Date:               16.12.2012
# Last Modified Date: 04.01.2022
# Last Modified By:   jn <jnordin@physik.hu-berlin.de>

import logging
from typing import Any, Union, Generator, Sequence, Literal, Optional
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib import cm
import pandas as pd
import numpy as np
import seaborn as sns
from adjustText import adjust_text

from ampel.types import UBson, T3Send
from ampel.struct.UnitResult import UnitResult
from ampel.view.T3Store import T3Store
from ampel.view.TransientView import TransientView
from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit
from ampel.util.pretty import prettyjson
from ampel.ztf.util.ZTFIdMapper import to_ztf_id



class HealpixCorrPlotter(AbsPhotoT3Unit):
	"""
	Compare healpix coordinate P-value with output from T2RunSncosmo.
	"""

	sncosmo_unit: str = 'T2RunSncosmo'
	model_name: Optional[str]   # Only use this model
	time_parameter: str = 't0'  # Name of the model parameter determining explosion / peak time

	# What do we study
	target_property: Literal['Abs fit peak mag', r'$\chi^2$ / d.o.f.'] = 'Abs fit peak mag'
	target_range: list[float] = [-13.5,-17.5]
	max_pvalue: float = 0.9

	# Plot params
	plotsize: Sequence[float] = [6,4]
    # List of inclusive lower limit, non-inc upper limit, marker type, label
#	marker_colors: list[str] = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]
#	marker_colors: list[str] = ["#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"]
	marker_colors = cm.get_cmap('plasma', 5).colors
#	marker_colors = [cm.get_cmap('summer', 3)(i) for i in range(3)]
	background_color: str = "tab:green"
	ndof_marker: list[Any] = [ [0,0.5,'o', marker_colors[0], '0 dof'], [1,1.5,'^', marker_colors[1], '1 dof'], [2,np.inf,'s', marker_colors[2], '>1 dof'], ]



	def process(self, gen: Generator[TransientView, T3Send, None], t3s: Optional[T3Store] = None) -> Union[UBson, UnitResult]:

		self.logger.info("Printing transients info")
		self.logger.info("=" * 80)
		count = 0

		table_rows: list[dict[str, Any]] = []
		for tran_view in gen:
			count += 1
			self.logger.info(str(count))
			# Stock info
			tinfo = self._get_stock_info(tran_view)

			# t2_info
			t2docs = tran_view.get_raw_t2_body(unit=self.sncosmo_unit)
			if t2docs is None:
				continue
			for t2info in t2docs:
				assert isinstance(t2info, dict)
				if self.model_name and not t2info['model_name']==self.model_name:
					continue
				tinfo['z'] = t2info['z']
				if t2info['z_source'] in ['AMPELz_group0', 'AMPELz_group1','AMPELz_group2','AMPELz_group3']:
					tinfo['z_sharp'] = True
				else:
					tinfo['z_sharp'] = False
				tinfo['zsource'] = t2info['z']
				tinfo['model'] = t2info['model_name']
				tinfo['model_peak_abs'] = t2info['fit_metrics']['restpeak_model_absmag_B']
				tinfo['model_peak_obs'] = t2info['fit_metrics']['obspeak_model_B']
				tinfo['ndof'] = t2info['sncosmo_result']['ndof']
				tinfo['chisq'] = t2info['sncosmo_result']['chisq']
				if t2info['sncosmo_result']['ndof']>0:
					tinfo['chisqndof'] = t2info['sncosmo_result']['chisq'] / t2info['sncosmo_result']['ndof']
				else:
					tinfo['chisqndof'] = -1.
				tinfo['time'] = t2info['sncosmo_result']['paramdict'][self.time_parameter]
			self.logger.info(tinfo)
			table_rows.append(tinfo)
		self.logger.info("=" * 80)
		self.logger.info(f"Printed info for {count} transients")


		df = pd.DataFrame.from_dict(table_rows)

		# figure
		if self.target_property=='Abs fit peak mag':
			df['target'] = df['model_peak_abs']
		elif self.target_property==r'$\chi^2$ / d.o.f.':
			df['target'] = df['chisqndof']
		dy = self.target_range[1] - self.target_range[0]

		fig = plt.figure(figsize=self.plotsize, dpi=300)
		ax = plt.gca()



		plt.fill_between([0,self.max_pvalue], self.target_range[0], self.target_range[1], alpha=0.5, color=self.background_color)
		plt.ylim([self.target_range[0]-dy, self.target_range[1]+dy])

        # Iterate through all detections. Cumbersome, but could find no way to make batch Plotting
        # while retaining PathCollection for adjustText and changing marker type.
		annotations = []
		targetpoints = []
		for i, row in df.iterrows():
            # Determine marker type (decided by d.o.f.)
			marker = None
			for minfo in self.ndof_marker:
				if (minfo[0] <= row['ndof'] < minfo[1]):
					marker = minfo[2]
					color = minfo[3]
			if marker is None:
			    continue
            # Determine outline (decided by redshift origin)
			markeredgecolor='None'
			if row['z_sharp']==True:
				markeredgecolor='k'
            # Determine whether to annotate (decided target region)
			if ((row['pvalue']<=self.max_pvalue) & (np.abs(row['target'])>=np.abs(self.target_range[0])) & (np.abs(row['target'])<=np.abs(self.target_range[1])) ):
				annotations.append( plt.text(row['pvalue'], row['target'], row['name']) )
				targetpoints.append( plt.plot(row['pvalue'], row['target'], marker, ms=10, color=color, markeredgecolor=markeredgecolor))
			else:
				plt.plot(row['pvalue'], row['target'], marker, ms=10, color=color, markeredgecolor=markeredgecolor, alpha=0.3)



		plt.xscale('log')
		ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
		ax.xaxis.set_minor_formatter(ticker.ScalarFormatter())
#		adjust_text(annotations, arrowprops=dict(arrowstyle='->', color='red'), add_objects=targetpoints)
		adjust_text(annotations, arrowprops=dict(arrowstyle='->', color='darkslategrey'))


        # Create legend
        #legend_elements = []
        #for minfo in self.ndof_marker:
        #    legend_elements.append( Line2D([0], [0], marker=minfo[2], color=minfo[3], label=minfo[4], markeredgecolor='None' ) )
		legend_elements = [Line2D([0], [0], marker=minfo[2], color=minfo[3], label=minfo[4], markeredgecolor='None', linewidth=0 ) for minfo in self.ndof_marker ]
		legend_elements.append(
            Line2D([0], [0], marker='o', markerfacecolor='None', label='Good z', markeredgecolor='k', linewidth=0 )
            )
		legend_elements.append(
            Patch(facecolor=self.background_color, edgecolor='None', label='Target region')
            )
		ax.legend(handles=legend_elements, loc='best', ncol=2)
#		plt.legend(loc=3,ncol=2)
		plt.xlabel('Healpix spatial P-value')
		plt.ylabel(self.target_property)
		# Determine titel (could be multiple?) from lists
		channels = set()
		for chlist in list(df['channel']):
			for ch in chlist:
				channels.add(ch)
		plt.title('{}'.format(' '.join(channels)))

		plt.savefig('/home/jnordin/tmp/test.pdf')

		return None


	def _get_stock_info(self, tran: TransientView) -> dict[Any,Any]:
		"""
		Gather relevant information from stock document.
		"""

		assert isinstance(tran.id, int)
		assert tran.stock
		stockinfo = {'id':tran.id, 'name': to_ztf_id(tran.id), 'channel': tran.stock['channel']}

		# Could there be multiple healpix journal entries? I guess it cannot be ruled out
		# FIXME: extra info should probably be in .extra, not mixed into the top level of JournalRecord
		hpixs = [el['healpix'] for el in tran.stock['journal'] if 'healpix' in el.keys()]
		if len(hpixs)==0:
			self.logger.info('No healpix info')
			return stockinfo
		stockinfo.update( hpixs[0])
		if len(hpixs)>1:
			for i, hpix in enumerate(hpixs[1:]):
				stockinfo.update({k+str(i):v for k,v in hpix.items()})
			self.logger.debug('Multiple healpix matches', extra=stockinfo)

		return stockinfo
