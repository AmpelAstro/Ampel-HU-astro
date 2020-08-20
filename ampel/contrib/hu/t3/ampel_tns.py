#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t3/ampel_tns.py
# License           : BSD-3-Clause
# Author            : Ken Smith
# Date              : May 2016
# Last Modified Date: Feb 2018
# Last Modified By  : Jakob Nordin

# -----------------------------------------------------------------------------
# A python sample code for sending a bulk report to the TNS.
# Original sample code by Ken Smith (May 2016) modified by Jakob Nordin (Feb 2018)
# -----------------------------------------------------------------------------

import re, requests, json, time
from collections import OrderedDict

TNSFILTERID = {1: "110", 2: "111", 3: "112"}
AT_REPORT_FORM = "bulk-report"
AT_REPORT_REPLY = "bulk-report-reply"
TNS_ARCHIVE = {'OTHER': '0', 'SDSS': '1', 'DSS': '2'}
TNS_BASE_URL_SANDBOX = "https://sandbox-tns.weizmann.ac.il/api/"
TNS_BASE_URL_REAL = "https://wis-tns.weizmann.ac.il/api/"

httpErrors = {
	304: 'Error 304: Not Modified: There was no new data to return.',
	400: 'Error 400: Bad Request: The request was invalid. An accompanying error message will explain why.',
	403: 'Error 403: Forbidden: The request is understood, but it has been refused. An accompanying error message will explain why',
	404: 'Error 404: Not Found: The URI requested is invalid or the resource requested, such as a category, does not exists.',
	500: 'Error 500: Internal Server Error: Something is broken.',
	503: 'Error 503: Service Unavailable.'
}


class TNSClient:
	"""Send Bulk TNS Request."""

	def __init__(self, baseURL, logger, options = {}):
		"""
		:param baseURL: Base URL of the TNS API
		:param options: (Default value = {})
		"""

		#self.baseAPIUrl = TNS_BASE_URL_SANDBOX
		self.baseAPIUrl = baseURL
		self.generalOptions = options
		self.logger = logger


	def buildUrl(self, resource):
		"""
		Build the full URL

		:param resource: the resource requested
		:return complete URL
		"""
		return self.baseAPIUrl + resource


	def buildParameters(self, parameters = {}):
		"""
		Merge the input parameters with the default parameters created when
		the class is constructed.

		:param parameters: input dict (Default value = {})
		:return p: merged dict
		"""
		p = self.generalOptions.copy()
		p.update(parameters)
		return p


	def jsonResponse(self, r):
		"""
		Send JSON response given requests object. Should be a python dict.

		:param r: requests object - the response we got back from the server
		:return d: json response converted to python dict
		"""

		d = {}
		# What response did we get?
		message = None
		status = r.status_code

		if status != 200:
			try:
				message = httpErrors[status]
			except ValueError:
				message = f'Error {status}: Undocumented error'

		if message is not None:
			self.logger.warn('TNS bulk submit: ' + message)
			return d

		# Did we get a JSON object?
		try:
			d = r.json()
		except ValueError as e:
			self.logger.error('TNS bulk submit: ' + e)
			d = {}
			return d

		# If so, what error messages if any did we get?
		self.logger.info(json.dumps(d, indent=4, sort_keys=True))

		if 'id_code' in d.keys() and 'id_message' in d.keys() and d['id_code'] != 200:
			self.logger.info(
				"TNS bulk submit: Bad response: code = %d, error = '%s'" %
				(d['id_code'], d['id_message'])
			)
		return d


	def sendBulkReport(self, options) -> dict:
		"""
		Send the JSON TNS request
		:param options: the JSON TNS request
		"""
		feed_url = self.buildUrl(AT_REPORT_FORM)
		feed_parameters = self.buildParameters({'data': (None, json.dumps(options))})

		# The requests.post method needs to receive a "files" entry, not "data".  And the "files"
		# entry needs to be a dictionary of tuples.  The first value of the tuple is None.
		self.logger.info('TNS bulk submit: ' + 'sending request')
		r = requests.post(feed_url, files = feed_parameters, timeout = 300)
		# Construct the JSON response and return it.
		self.logger.info('TNS bulk submit: ' + 'got response (or timed out)')
		return self.jsonResponse(r)


	def bulkReportReply(self, options):
		"""
		Get the report back from the TNS

		:param options: dict containing the report ID
		:return: dict

		"""
		feed_url = self.buildUrl(AT_REPORT_REPLY)
		feed_parameters = self.buildParameters(options)

		self.logger.info('TNS bulk submit: ' + 'looking for reply report')
		r = requests.post(feed_url, files = feed_parameters, timeout = 300)
		self.logger.info('TNS bulk submit: ' + 'got report (or timed out)')
		return self.jsonResponse(r)


def addBulkReport(report, tnsApiKey, logger, sandbox=True):
	"""
	Send the report to the TNS

	:param tnsBaseURL: TNS base URL
	:param tnsApiKey: TNS API Key
	:param sandbox: Submits to TNS sandbox area if true
	:return reportId: TNS report ID

	"""

	if sandbox:
		baseurl = TNS_BASE_URL_SANDBOX
		logger.warning("TNS report submitted to sandbox")
	else:
		baseurl = TNS_BASE_URL_REAL

	feed_handler = TNSClient(baseurl, logger, {'api_key': (None, tnsApiKey)})
	reply = feed_handler.sendBulkReport(report)

	reportId = None

	if reply:
		try:
			reportId = reply['data']['report_id']
			logger.info('TNS bulk submit: successful with ID %s' % (reportId))
		except ValueError:
			logger.error("Empty response. Something went wrong. Is the API Key OK?")
		except KeyError:
			logger.error("Cannot find the data key. Something is wrong.")

	return reportId


def getBulkReportReply(reportId, tnsApiKey, logger, sandbox=True):
	"""
	Get the TNS response for the specified report ID
	:param tnsApiKey: TNS API Key
	:return request: The original request
	:return response: The TNS response
	"""

	if sandbox:
		baseurl = TNS_BASE_URL_SANDBOX
		logger.warning("TNS report submitted to sandbox")
	else:
		baseurl = TNS_BASE_URL_REAL

	feed_handler = TNSClient(baseurl, logger, {'api_key': (None, tnsApiKey)})
	reply = feed_handler.bulkReportReply({'report_id': (None, str(reportId))})

	request = None
	response = None
	# reply should be a dict
	if (reply and 'id_code' in reply.keys() and reply['id_code'] == 404):
		logger.warn(
			f"TNS bulk submit {reportId}: Unknown report. "
			f"Perhaps the report has not yet been processed."
		)


	if reply and 'id_code' in reply.keys() and reply['id_code'] == 200:
		try:
			request = reply['data']['received_data']['at_report']
			response = reply['data']['feedback']['at_report']
		except ValueError:
			logger.error("TNS bulk submit: cannot find the response feedback payload.")

	logger.info(f'TNS bulk submit: got response {response} for request {request}')

	return response


# function for search obj
def tnsName(ra, dec, api_key, sandbox=True, matchradius=5):
	"""
	Search for TNS transients around input ra, dec (given in degrees)
	Returns None or interl name
	"""

	if sandbox:
		baseurl = TNS_BASE_URL_SANDBOX
	else:
		baseurl = TNS_BASE_URL_REAL
	search_url = baseurl + 'get/search'

	search_obj = {"ra": ra, "dec": dec, "radius": matchradius, "units": "arcsec"}
	search_data = [
		('api_key', (None, api_key)),
		('data', (None, json.dumps(search_obj)))
	]

	# Do request
	response = requests.post(search_url, files=search_data)

	# Verify connection
	status = response.status_code
	if status != 200:
		try:
			message = httpErrors[status]
		except ValueError:
			message = f'Error {status}: Undocumented error'
		return False, message

	# Did we get a reply?
	if response is None:
		return False, 'Error: No TNS response'

	# Can we find a single name?
	parsed = json.loads(response.text, object_pairs_hook=OrderedDict)
	try:
		tnsnames = [v['prefix'] + v['objname'] for v in parsed["data"]["reply"]]
	except KeyError:
		return False, 'Error: No TNS names in response'

	return tnsnames, 'Found TNS name(s)'


def tnsInternal(tnsname, api_key, sandbox=True):
	"""
	Get internal names for a given TNS transient.
	Required to avoid oversubmitting transients?
	"""

	if sandbox:
		baseurl = TNS_BASE_URL_SANDBOX
	else:
		baseurl = TNS_BASE_URL_REAL
	get_url = baseurl + 'get/object'

	get_obj = OrderedDict([("objname", tnsname), ("photometry", "0"), ("spectra", "0")])
	get_data = [('api_key', (None, api_key)), ('data', (None, json.dumps(get_obj)))]

	# get obj using request module
	response = requests.post(get_url, files=get_data)

	# Verify connection
	status = response.status_code
	if status != 200:
		try:
			message = httpErrors[status]
		except ValueError:
			message = f'Error {status}: Undocumented error'
		return False, message

	# Did we get a reply?
	if response is None:
		return False, 'Error: No TNS response'

	parsed = json.loads(response.text, object_pairs_hook=OrderedDict)
	if parsed['data']['reply']['internal_name'] is None:
		return None, 'No internal TNS name'

	return parsed['data']['reply']['internal_name'], 'Got internal name response'


# function for search obj
def tnssearch(json_list, api_key, sandbox=True):

	if sandbox:
		baseurl = TNS_BASE_URL_SANDBOX
		#baseurl="https://sandbox-tns.weizmann.ac.il/api/get"
	else:
		baseurl = TNS_BASE_URL_REAL
		#baseurl="https://wis-tns.weizmann.ac.il/api/get"

	try:
		# url for search obj
		search_url = baseurl + 'get/search'
		# change json_list to json format
		json_file = OrderedDict(json_list)
		# construct the list of (key, value) pairs
		search_data = [('api_key', (None, api_key)), ('data', (None, json.dumps(json_file)))]
		# search obj using request module
		response = requests.post(search_url, files=search_data)
		return response
	except Exception as e:
		return [None, 'Error message : \n' + str(e)]


# function for get obj
def tnsget(json_list, api_key, sandbox=True):

	if sandbox:
		baseurl = TNS_BASE_URL_SANDBOX
		#baseurl="https://sandbox-tns.weizmann.ac.il/api/get"
	else:
		baseurl = TNS_BASE_URL_REAL
		#baseurl="https://wis-tns.weizmann.ac.il/api/get"

	try:
		# url for get obj
		get_url = baseurl + 'get/object'
		# change json_list to json format
		json_file = OrderedDict(json_list)
		# construct the list of (key, value) pairs
		get_data = [('api_key', (None, api_key)), ('data', (None, json_file))]
		# get obj using request module
		response = requests.post(get_url, files=get_data)
		return response
	except Exception as e:
		return [None, 'Error message : \n' + str(e)]


def sendTNSreports(atreportlists, api_key, logger, sandbox=True):
	"""
	Based on a lists of reportlists, send to TNS.
	Return results for journal entries
	"""

	# Submit to TNS
	MAX_LOOP = 25
	SLEEP = 2

	reportresult = {}
	for atreport in atreportlists:

		# Submit a report
		for _ in range(MAX_LOOP):
			reportid = addBulkReport(atreport, api_key, logger, sandbox=sandbox)
			if reportid:
				logger.info('TNS report ID %s' % (reportid))
				break
			time.sleep(SLEEP)
		else:
			logger.info('TNS bulk report failed')
			continue

		# Try to read reply
		for _ in range(MAX_LOOP):
			time.sleep(SLEEP)
			response = getBulkReportReply(reportid, api_key, logger, sandbox=sandbox)
			if isinstance(response, list):
				break
		else:
			logger.info("TNS Report reading failed")
			continue

		# Parse reply for evaluation
		for k, v in atreport["at_report"].items():
			if '100' in response[k].keys():
				logger.info("TNS Inserted with name %s" % (response[k]["100"]["objname"]))
				reportresult[v["internal_name"]] = ['TNS inserted', {"TNSName": response[k]["100"]["objname"]}]
			elif '101' in response[k].keys():
				logger.info("ALready existing with name %s" % (response[k]["101"]["objname"]))
				reportresult[v["internal_name"]] = ['TNS pre-existing', {"TNSName": response[k]["101"]["objname"]}]

	return reportresult


def get_tnsname(ra, dec, api_key, logger, sandbox=True):
	"""
	look for names registered at tns for a given position
	"""
	tns_name = None

	# Look for TNS name at the coordinate of the transient
	tnsnames, runstatus = tnsName(ra, dec, api_key, sandbox=sandbox)
	if re.match('Error', runstatus):
		logger.info("TNS get error", extra={"tns_request": runstatus})
		return None, []
	if len(tnsnames) > 1:
		logger.debug("Multipe TNS names, choosing first", extra={"tns_names": tnsnames})
		tns_name = tnsnames[0]
	elif len(tnsnames) == 1:
		tns_name = tnsnames[0]
	elif len(tnsnames) == 0:
		# No TNS name, then no need to look for internals
		return tns_name, None
	logger.info("TNS get cand id", extra={"tns_name": tns_name})

	# Look for internal name (note that we skip the prefix)
	internal_name, runstatus = tnsInternal(tns_name[2:], api_key, sandbox=sandbox)
#		if internal_name is not None:
#			self.logger.info("TNS search", extra={"tns_internal_name": internal_name})
#			pass
#			if re.search('ZTF', internal_name):
#				if not internal_name == sne[1]["ztf_name"]:
#					#self.logger.info("TNS registered under other ZTF name %s", extra={"tns_other_internal": internal_name})
	return tns_name, internal_name
