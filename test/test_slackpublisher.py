
import pytest
from ampel.contrib.hu.t3.SlackSummaryPublisher import SlackSummaryPublisher
pytest_plugins = ['ampel.test.fixtures']

import requests, slackclient
import csv
from io import StringIO

def test_slacksummary(t3_transient_views, mocker):

	run_config = {
		"mycols": [
			"ztf_name",
			"ra",
			"dec",
			"magpsf",
			"sgscore1",
			"rb",
			"most_recent_detection",
			"first_detection",
			"n_detections",
			"distnr",
			"distpsnr1",
			"isdiffpos",
			"_id"
		],
		"channel(s)": [ '0', '1'],
		"excitement_levels": {
			"Low": 50,
			"Mid": 200,
			"High": 400
		},
		"Slack_token": "xoxoxox",
		"Slack_channel": "#ampel-live",
		"full_photometry": True
	}
	
	assert len(t3_transient_views) < run_config['excitement_levels']['Low'], 'Small number passed'

	t3 = SlackSummaryPublisher(None, run_config=run_config)

	# intercept Slack API calls
	mocker.patch('requests.post')
	mocker.patch('slackclient.SlackClient.api_call')

	t3.add(t3_transient_views)
	t3.done()
	
	api_call = slackclient.SlackClient.api_call
	api_call.assert_called_once()
	assert 'MEH!' in api_call.call_args[1]['text'], 'Text matches number of transients selected'

	requests.post.assert_called()
	assert len(requests.post.call_args_list) == 2, '2 explicit requests issued'

	# Verify summary
	content = requests.post.call_args_list[0][1]['files']['file']
	with StringIO(content) as f:
		reader = csv.DictReader(f)
		rows = list(reader)

	# verify that T2 information is in summary
	t2s = set(t2.t2_unit_id for tv in t3_transient_views for t2 in tv.t2records)
	assert len(t2s) > 0

        # Verify that nested t2 results were extrected
        # This tests assumes that sncosmo was run on the test data
	for key in ["T2-model","T2-sncosmo_info_ndof","T2-fit_results_t0"]:
		assert key in reader.fieldnames

	assert len(rows) == len(t3_transient_views), '1 row per transient'

	# Verify photometry dump
	content = requests.post.call_args_list[1][1]['files']['file']
	with StringIO(content) as f:
		rows = list(csv.reader(f))
	assert len(rows) == sum(len(v.photopoints) for v in t3_transient_views)+1, '1 row per photopoint'


