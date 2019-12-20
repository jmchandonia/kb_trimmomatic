# -*- coding: utf-8 -*-
#BEGIN_HEADER
import sys
import traceback
from biokbase.workspace.client import Workspace as workspaceService
import requests
requests.packages.urllib3.disable_warnings()
import subprocess
import os
import re
from pprint import pprint, pformat
from datetime import datetime
import uuid

## SDK Utils
from ReadsUtils.ReadsUtilsClient import ReadsUtils
from SetAPI.SetAPIServiceClient import SetAPI
from DataFileUtil.DataFileUtilClient import DataFileUtil as DFUClient
from KBaseReport.KBaseReportClient import KBaseReport
#END_HEADER


class kb_trimmomatic:
    '''
    Module Name:
    kb_trimmomatic

    Module Description:
    A KBase module: kb_trimmomatic
This module contains two methods

runTrimmomatic() to backend a KBase App, potentially operating on ReadSets
execTrimmomatic() the local method that handles overloading Trimmomatic to run on a set or a single library
execTrimmomaticSingleLibrary() runs Trimmomatic on a single library
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "1.2.11"
    GIT_URL = "https://github.com/kbaseapps/kb_trimmomatic"
    GIT_COMMIT_HASH = "9ff31d23f62491d7a47c004f6cf8f800535b47f1"

    #BEGIN_CLASS_HEADER
    workspaceURL = None
    TRIMMOMATIC = 'java -jar /kb/module/Trimmomatic-0.36/trimmomatic-0.36.jar'
    ADAPTER_DIR = '/kb/module/Trimmomatic-0.36/adapters/'

    def log(self, target, message):
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()


    # Determine if phred64
    #
    def is_fastq_phred64 (self, this_input_path):
        read_buf_size  = 65536
        input_is_phred33 = False
        data_seen = False
        with open (this_input_path, 'r', read_buf_size) as this_input_handle:
            while True:
                line = this_input_handle.readline()
                if not line:
                    break
                if not line.startswith('@'):
                    raise ValueError ("Badly formatted FASTQ file: "+this_input_path+"\n"+"BAD LINE: '"+line+"'")
                # skip two more lines
                this_input_handle.readline()  # seq
                this_input_handle.readline()  # '+' qual header

                qual_line = this_input_handle.readline().rstrip()
                data_seen = True
                #def qual33(qual64): return chr(ord(qual64)-31)
                for qual_val in qual_line:
                    q64_ascii = ord(qual_val)
                    if q64_ascii < 64:
                        input_is_phred33 = True
                        break
                if input_is_phred33:
                    break
        if not data_seen:
            raise ValueError ("no qual score line found in FASTQ file: "+this_input_path)

        input_is_phred64 = not input_is_phred33
        return input_is_phred64


    # Translate phred64 to phred33
    #
    def translate_fastq_from_phred64_to_phred33 (self, this_input_path, this_output_path):

        if not self.is_fastq_phred64 (this_input_path):
            return this_input_path

        # internal Method
        def qual33(qual64): return chr(ord(qual64)-31)

        # read through and translate qual scores
        read_buf_size  = 65536
        write_buf_size = 65536

        qual33_handle = open (this_output_path, 'w', write_buf_size)
        with open (this_input_path, 'r', read_buf_size) as this_input_handle:
            while True:
                buf = []
                line = this_input_handle.readline()
                if not line:
                    break
                if line.startswith('@'):
                    buf.append(line)  # header
                    buf.append(this_input_handle.readline())  # seq
                    buf.append(this_input_handle.readline())  # '+'

                    qual_line = this_input_handle.readline().rstrip()
                    q33_line = ''
                    for q64 in qual_line:
                        q33_line += qual33(q64)
                    buf.append(q33_line+"\n")
                    qual33_handle.write(''.join(buf))
        qual33_handle.close()

        return this_output_path


    # Set up Trimmomatic params
    #
    def parse_trimmomatic_steps(self, input_params):
        # validate input parameters and return string defining trimmomatic steps

        parameter_string = ''

#        if 'read_type' not in input_params and input_params['read_type'] is not None:
#            raise ValueError('read_type not defined')
#        elif input_params['read_type'] not in ('PE', 'SE'):
#            raise ValueError('read_type must be PE or SE')

#        if 'quality_encoding' not in input_params and input_params['quality_encoding'] is not None:
#            raise ValueError('quality_encoding not defined')
#        elif input_params['quality_encoding'] not in ('phred33', 'phred64'):
#            raise ValueError('quality_encoding must be phred33 or phred64')

        # set adapter trimming
        if ('adapterFa' in input_params and input_params['adapterFa'] is not None and
            'seed_mismatches' in input_params and input_params['seed_mismatches'] is not None and
            'palindrome_clip_threshold' in input_params and input_params['palindrome_clip_threshold'] is not None and
            'simple_clip_threshold' in input_params and input_params['simple_clip_threshold'] is not None):
            parameter_string = ("ILLUMINACLIP:" + self.ADAPTER_DIR +
                                ":".join((str(input_params['adapterFa']),
                                          str(input_params['seed_mismatches']),
                                          str(input_params['palindrome_clip_threshold']),
                                          str(input_params['simple_clip_threshold']))) + " ")
        elif ( ('adapterFa' in input_params and input_params['adapterFa'] is not None) or
               ('seed_mismatches' in input_params and input_params['seed_mismatches'] is not None) or
               ('palindrome_clip_threshold' in input_params and input_params['palindrome_clip_threshold'] is not None) or
               ('simple_clip_threshold' in input_params and input_params['simple_clip_threshold'] is not None) ):
            raise ValueError('Adapter Clipping requires Adapter, Seed Mismatches, Palindrome Clip Threshold and Simple Clip Threshold')

        # set Crop
        if 'crop_length' in input_params and input_params['crop_length'] is not None \
                and int(input_params['crop_length']) > 0:
            parameter_string += 'CROP:' + str(input_params['crop_length']) + ' '

        # set Headcrop
        if 'head_crop_length' in input_params and input_params['head_crop_length'] is not None \
                and input_params['head_crop_length'] > 0:
            parameter_string += 'HEADCROP:' + str(input_params['head_crop_length']) + ' '

        # set Leading
        if 'leading_min_quality' in input_params and input_params['leading_min_quality'] is not None \
                and input_params['leading_min_quality'] > 0:
            parameter_string += 'LEADING:' + str(input_params['leading_min_quality']) + ' '

        # set Trailing
        if 'trailing_min_quality' in input_params and input_params['trailing_min_quality'] is not None \
                and input_params['trailing_min_quality'] > 0:
            parameter_string += 'TRAILING:' + str(input_params['trailing_min_quality']) + ' '

        # set sliding window
        if 'sliding_window_size' in input_params and input_params['sliding_window_size'] is not None \
                and input_params['sliding_window_size'] > 0 \
                and 'sliding_window_min_quality' in input_params and input_params['sliding_window_min_quality'] is not None \
                and input_params['sliding_window_min_quality'] > 0:
            parameter_string += 'SLIDINGWINDOW:' + str(input_params['sliding_window_size']) + ":" + str(input_params['sliding_window_min_quality']) + ' '
        elif ('sliding_window_size' in input_params and input_params['sliding_window_size'] is not None \
                and input_params['sliding_window_size'] > 0) \
             or ('sliding_window_min_quality' in input_params and input_params['sliding_window_min_quality'] is not None \
                and input_params['sliding_window_size'] > 0):
            raise ValueError('Sliding Window filtering requires both Window Size and Window Minimum Quality to be set')

        # set min length
        if 'min_length' in input_params and input_params['min_length'] is not None \
                and input_params['min_length'] > 0:
            parameter_string += 'MINLEN:' + str(input_params['min_length']) + ' '

        if parameter_string == '':
            raise ValueError('No filtering/trimming steps specified!')

        return parameter_string

    def _save_RNASeqSampleSet(self, items, wsName, output_SampleSet_name, reads_desc_ext,
                              single_reads):

        print ('Start saving RNASeqSampleSet object')
        workspace_id = self.dfu.ws_name_to_id(wsName)
        Library_type = 'SingleEnd' if single_reads else 'PairedEnd'
        sample_set_data = {'sampleset_id': output_SampleSet_name,
                           'sampleset_desc': reads_desc_ext,
                           'Library_type': Library_type,
                           'sample_ids': [item.get('ref') for item in items],
                           'condition': [item.get('label') for item in items],
                           'domain': 'Unknown',
                           'num_samples': len(items),
                           'platform': 'Unknown'}
        save_object_params = {
            'id': workspace_id,
            'objects': [{'type': 'KBaseRNASeq.RNASeqSampleSet',
                         'data': sample_set_data,
                         'name': output_SampleSet_name}]
        }

        dfu_oi = self.dfu.save_objects(save_object_params)[0]
        sample_set_ref = str(dfu_oi[6]) + '/' + str(dfu_oi[0]) + '/' + str(dfu_oi[4])

        return sample_set_ref

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.scratch = os.path.abspath(config['scratch'])
        self.handleURL = config['handle-service-url']
        self.serviceWizardURL = config['service-wizard-url']

        self.callbackURL = os.environ.get('SDK_CALLBACK_URL', None)
        if self.callbackURL is None:
            raise ValueError("SDK_CALLBACK_URL not set in environment")

        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)
        os.chdir(self.scratch)

        self.dfu = DFUClient(self.callbackURL)
        #END_CONSTRUCTOR
        pass


    def runTrimmomatic(self, ctx, input_params):
        """
        :param input_params: instance of type "runTrimmomaticInput"
           (runTrimmomatic() ** ** to backend a KBase App, potentially
           operating on ReadSets) -> structure: parameter "input_ws" of type
           "workspace_name" (** Common types), parameter "input_reads_ref" of
           type "data_obj_ref", parameter "output_ws" of type
           "workspace_name" (** Common types), parameter "output_reads_name"
           of type "data_obj_name", parameter "translate_to_phred33" of type
           "bool", parameter "adapter_clip" of type "AdapterClip_Options" ->
           structure: parameter "adapterFa" of String, parameter
           "seed_mismatches" of Long, parameter "palindrome_clip_threshold"
           of Long, parameter "simple_clip_threshold" of Long, parameter
           "sliding_window" of type "SlidingWindow_Options" (parameter
           groups) -> structure: parameter "sliding_window_size" of Long,
           parameter "sliding_window_min_quality" of Long, parameter
           "leading_min_quality" of Long, parameter "trailing_min_quality" of
           Long, parameter "crop_length" of Long, parameter
           "head_crop_length" of Long, parameter "min_length" of Long
        :returns: instance of type "runTrimmomaticOutput" -> structure:
           parameter "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN runTrimmomatic
        console = []
        self.log(console, 'Running runTrimmomatic with parameters: ')
        self.log(console, "\n"+pformat(input_params))

        token = ctx['token']
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token

        SERVICE_VER = 'release'

        # param checks
        if ('output_ws' not in input_params or input_params['output_ws'] is None):
            input_params['output_ws'] = input_params['input_ws']

        required_params = ['input_reads_ref',
                           'output_ws',
                           'output_reads_name'
#                           'read_type'
                          ]
        for required_param in required_params:
            if required_param not in input_params or input_params[required_param] == None:
                raise ValueError ("Must define required param: '"+required_param+"'")

        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        provenance[0]['input_ws_objects']=[str(input_params['input_reads_ref'])]

        # set up and run execTrimmomatic()
        #
        execTrimmomaticParams = { 'input_reads_ref': str(input_params['input_reads_ref']),
                                  'output_ws': input_params['output_ws'],
                                  'output_reads_name': input_params['output_reads_name']
#                                  'read_type': input_params['read_type'],
                                 }

        #if 'quality_encoding' in input_params:
        #    execTrimmomaticParams['quality_encoding'] = input_params['quality_encoding']
        if 'translate_to_phred33' in input_params:
            execTrimmomaticParams['translate_to_phred33'] = input_params['translate_to_phred33']

        # adapter_clip grouped params
        if 'adapter_clip' in input_params and input_params['adapter_clip'] != None:
            if 'adapterFa' in input_params['adapter_clip']:
                execTrimmomaticParams['adapterFa'] = input_params['adapter_clip']['adapterFa']
            else:
                execTrimmomaticParams['adapterFa'] = None

            if 'seed_mismatches' in input_params['adapter_clip']:
                execTrimmomaticParams['seed_mismatches'] = input_params['adapter_clip']['seed_mismatches']
            else:
                execTrimmomaticParams['seed_mismatches'] = None

            if 'palindrome_clip_threshold' in input_params['adapter_clip']:
                execTrimmomaticParams['palindrome_clip_threshold'] = input_params['adapter_clip']['palindrome_clip_threshold']
            else:
                execTrimmomaticParams['palindrome_clip_threshold'] = None

            if 'simple_clip_threshold' in input_params['adapter_clip']:
                execTrimmomaticParams['simple_clip_threshold'] = input_params['adapter_clip']['simple_clip_threshold']
            else:
                execTrimmomaticParams['simple_clip_threshold'] = None

        # sliding window
        if 'sliding_window' in input_params:
            if 'sliding_window_size' in input_params['sliding_window']:
                execTrimmomaticParams['sliding_window_size'] = input_params['sliding_window']['sliding_window_size']
            else:
                execTrimmomaticParams['sliding_window_size'] = None

            if 'sliding_window_min_quality' in input_params['sliding_window']:
                execTrimmomaticParams['sliding_window_min_quality'] = input_params['sliding_window']['sliding_window_min_quality']
            else:
                execTrimmomaticParams['sliding_window_min_quality'] = None

        # remaining params
        if 'leading_min_quality' in input_params:
            execTrimmomaticParams['leading_min_quality'] = input_params['leading_min_quality']
        if 'trailing_min_quality' in input_params:
            execTrimmomaticParams['trailing_min_quality'] = input_params['trailing_min_quality']
        if 'crop_length' in input_params:
            execTrimmomaticParams['crop_length'] = input_params['crop_length']
        if 'head_crop_length' in input_params:
            execTrimmomaticParams['head_crop_length'] = input_params['head_crop_length']
        if 'min_length' in input_params:
            execTrimmomaticParams['min_length'] = input_params['min_length']

        # RUN
        trimmomatic_retVal = self.execTrimmomatic (ctx, execTrimmomaticParams)[0]


        # build report
        #
        reportName = 'kb_trimmomatic_report_'+str(uuid.uuid4())

        reportObj = {'objects_created': [],
                     #'text_message': '',  # or is it 'message'?
                     'message': '',  # or is it 'text_message'?
                     'direct_html': None,
                     'file_links': [],
                     'html_links': [],
                     'html_window_height': 220,
                     'workspace_name': input_params['input_ws'],
                     'report_object_name': reportName
                     }


        # text report (replaced by HTML report)
        try:
            #reportObj['text_message'] = trimmomatic_retVal['report']
            #reportObj['message'] = trimmomatic_retVal['report']
            msg = trimmomatic_retVal['report']
        except:
            raise ValueError ("no report generated by execTrimmomatic()")

        # parse text report
        report_data = []
        report_field_order = []
        report_lib_refs = []
        report_lib_names = []
        lib_i = -1

        # This is some powerful brute force nonsense, but it should be okay.
        #   (Note: it was not OK.  Now it is)
        se_expected_field_order = ['Input Reads',
                                   'Surviving',
                                   'Dropped']
        se_report_re = re.compile('^Input Reads:\s*(\d+)\s*Surviving:\s*(\d+)\s*\(\d+\.\d+\%\)\s*Dropped:\s*(\d+)\s*\(\d+\.\d+\%\)')
        pe_expected_field_order = ['Input Read Pairs',
                                   'Both Surviving',
                                   'Forward Only Surviving',
                                   'Reverse Only Surviving',
                                   'Dropped']
        pe_report_re = re.compile('^Input Read Pairs:\s*(\d+)\s*Both Surviving:\s*(\d+)\s*\(\d+\.\d+\%\)\s*Forward Only Surviving:\s*(\d+)\s*\(\d+\.\d+\%\)\s*Reverse Only Surviving:\s*(\d+)\s*\(\d+\.\d+\%\)\s*Dropped:\s*(\d+)\s*\(\d+\.\d+\%\)')

        for line in trimmomatic_retVal['report'].split("\n"):
            if line.startswith("RUNNING"):
                lib_i += 1
                lib_ids = re.sub("RUNNING TRIMMOMATIC ON LIBRARY: ", '', line)
                [ref, name] = lib_ids.split(" ")
                report_lib_refs.append(ref)
                report_lib_names.append(name)
                report_data.append({})
                report_field_order.append([])
            elif line.startswith("-"):
                continue
            elif len(line) == 0:
                continue
            else:
                # single end stats
                m_se = se_report_re.match(line)
                if m_se and len(m_se.groups()) == len(se_expected_field_order):
                    report_field_order[lib_i] = se_expected_field_order
                    report_data[lib_i] = dict(zip(report_field_order[lib_i], m_se.groups()))
                    for f_name in report_field_order[lib_i]:
                        report_data[lib_i][f_name] = int(report_data[lib_i][f_name])


                # paired end stats
                m_pe = pe_report_re.match(line)
                if m_pe and len(m_pe.groups()) == len(pe_expected_field_order):
                    report_field_order[lib_i] = pe_expected_field_order
                    report_data[lib_i] = dict(zip(report_field_order[lib_i], m_pe.groups()))
                    for f_name in report_field_order[lib_i]:
                        report_data[lib_i][f_name] = int(report_data[lib_i][f_name])


                else:
                    self.log(console, "SKIPPING OUTPUT.  Can't parse [" + line + "] (lib_i=" + str(lib_i) + ")")


        #### HTML report
        ##
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        html_output_dir = os.path.join(self.scratch,'output_html.'+str(timestamp))
        if not os.path.exists(html_output_dir):
            os.makedirs(html_output_dir)
        html_file = input_params['output_reads_name']+'.html'
        output_html_file_path = os.path.join(html_output_dir, html_file);

        # html config
        sp = '&nbsp;'
        text_color = "#606060"
        bar_color = "lightblue"
        bar_width = 100
        bar_char = "."
        bar_fontsize = "-2"
        row_spacing = "-2"

        html_report_lines = ['<html>']
        html_report_lines += ['<body bgcolor="white">']

        for lib_i in range(len(report_data)):
            html_report_lines += ['<p><b><font color="'+text_color+'">TRIMMOMATIC RESULTS FOR '+str(report_lib_names[lib_i])+' (object '+str(report_lib_refs[lib_i])+')</font></b><br>'+"\n"]
            high_val = 0
            if not len(report_field_order[lib_i]) > 0:
                html_report_lines += ['All reads were trimmed - no new reads object created.']
            else:
                html_report_lines += ['<table cellpadding=0 cellspacing=0 border=0>']
                html_report_lines += ['<tr><td></td><td>'+sp+sp+sp+sp+'</td><td></td><td>'+sp+sp+'</td></tr>']
                for f_name in report_field_order[lib_i]:
                    if int(report_data[lib_i][f_name]) > high_val:
                        high_val = int(report_data[lib_i][f_name])

                for f_name in report_field_order[lib_i]:
                    percent = round(float(report_data[lib_i][f_name])/float(high_val)*100, 1)
                    this_width = int(round(float(bar_width)*float(report_data[lib_i][f_name])/float(high_val), 0))

                    #self.log(console,"this_width: "+str(this_width)+" report_data: "+str(report_data[lib_i][f_name])+" calc: "+str(float(width)*float(report_data[lib_i][f_name])/float(high_val)))  # DEBUG
                    if this_width < 1:
                        if report_data[lib_i][f_name] > 0:
                            this_width = 1
                        else:
                            this_width = 0
                    html_report_lines += ['<tr>']
                    html_report_lines += ['    <td align=right><font color="'+text_color+'">'+str(f_name)+'</font></td><td></td>']
                    #html_report_lines += ['    <td align=right><font color="'+text_color+'">'+'{:0,}'.format(report_data[lib_i][f_name])+'</font></td><td></td>']
                    html_report_lines += ['    <td align=right><font color="'+text_color+'">'+str(report_data[lib_i][f_name])+'</font></td><td></td>']
                    html_report_lines += ['    <td align=right><font color="'+text_color+'">'+'('+str(percent)+'%)'+sp+sp+'</font></td><td></td>']

                    if this_width > 0:
                        for tic in range(this_width):
                            html_report_lines += ['    <td bgcolor="'+bar_color+'"><font size='+bar_fontsize+' color="'+bar_color+'">'+bar_char+'</font></td>']
                    html_report_lines += ['</tr>']
                    html_report_lines += ['<tr><td><font size='+row_spacing+'>'+sp+'</font></td></tr>']

                html_report_lines += ['</table>']
                html_report_lines += ['<p>']
        html_report_lines += ['</body>']
        html_report_lines += ['</html>']

        # write html to file and upload
        html_report_str = "\n".join(html_report_lines)
        #reportObj['direct_html'] = "\n".join(html_report_lines)   # doesn't always fit in buf
        with open (output_html_file_path, 'w', 0) as html_handle:
            html_handle.write(html_report_str)

        try:
            html_upload_ret = self.dfu.file_to_shock({'file_path': html_output_dir,
            #html_upload_ret = dfu.file_to_shock({'file_path': output_html_file_path,
                                                 #'make_handle': 0})
                                                 'make_handle': 0,
                                                 'pack': 'zip'})
        except:
            raise ValueError ('error uploading HTML file to shock')

        # attach to report obj
        #reportObj['direct_html'] = None
        reportObj['direct_html'] = ''
        reportObj['direct_html_link_index'] = 0
        reportObj['html_links'] = [{'shock_id': html_upload_ret['shock_id'],
                                    'name': html_file,
                                    'label': input_params['output_reads_name']+' HTML'
                                    }
                                   ]

        # trimmed object
        if trimmomatic_retVal['output_filtered_ref'] != None:
            try:
                # DEBUG
                #self.log(console,"OBJECT CREATED: '"+str(trimmomatic_retVal['output_filtered_ref'])+"'")

                reportObj['objects_created'].append({'ref':trimmomatic_retVal['output_filtered_ref'],
                                                     'description':'Trimmed Reads'})
            except:
                raise ValueError ("failure saving trimmed output")
        else:
            self.log(console, "No trimmed output generated by execTrimmomatic()")

        if trimmomatic_retVal.get('output_filtered_sampleset_ref'):
            try:
                reportObj['objects_created'].append({'ref':trimmomatic_retVal['output_filtered_sampleset_ref'],
                                                     'description':'Trimmed Reads'})
            except:
                raise ValueError ("failure saving trimmed output")
        else:
            self.log(console, "No trimmed output generated by execTrimmomatic()")

        # unpaired fwd
        if trimmomatic_retVal['output_unpaired_fwd_ref'] != None:
            try:
                reportObj['objects_created'].append({'ref':trimmomatic_retVal['output_unpaired_fwd_ref'],
                                                     'description':'Trimmed Unpaired Forward Reads'})
            except:
                raise ValueError ("failure saving unpaired fwd output")
        else:
            pass

        if trimmomatic_retVal.get('output_unpaired_sampleset_fwd_ref'):
            try:
                reportObj['objects_created'].append({'ref':trimmomatic_retVal['output_unpaired_sampleset_fwd_ref'],
                                                     'description':'Trimmed Unpaired Forward Reads'})
            except:
                raise ValueError ("failure saving unpaired fwd output")
        else:
            pass

        # unpaired rev
        if trimmomatic_retVal['output_unpaired_rev_ref'] != None:
            try:
                reportObj['objects_created'].append({'ref':trimmomatic_retVal['output_unpaired_rev_ref'],
                                                     'description':'Trimmed Unpaired Reverse Reads'})
            except:
                raise ValueError ("failure saving unpaired fwd output")
        else:
            pass

        if trimmomatic_retVal.get('output_unpaired_sampleset_rev_ref'):
            try:
                reportObj['objects_created'].append({'ref':trimmomatic_retVal['output_unpaired_sampleset_rev_ref'],
                                                     'description':'Trimmed Unpaired Reverse Reads'})
            except:
                raise ValueError ("failure saving unpaired fwd output")
        else:
            pass

        # save report object
        #
        report = KBaseReport(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
        #report_info = report.create({'report':reportObj, 'workspace_name':input_params['input_ws']})
        report_info = report.create_extended_report(reportObj)

        output = { 'report_name': report_info['name'], 'report_ref': report_info['ref'] }
        #END runTrimmomatic

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method runTrimmomatic return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def execTrimmomatic(self, ctx, input_params):
        """
        :param input_params: instance of type "execTrimmomaticInput"
           (execTrimmomatic() ** ** the local method that runs Trimmomatic on
           each read library) -> structure: parameter "input_reads_ref" of
           type "data_obj_ref", parameter "output_ws" of type
           "workspace_name" (** Common types), parameter "output_reads_name"
           of type "data_obj_name", parameter "read_type" of String,
           parameter "adapterFa" of String, parameter "seed_mismatches" of
           Long, parameter "palindrome_clip_threshold" of Long, parameter
           "simple_clip_threshold" of Long, parameter "translate_to_phred33"
           of type "bool", parameter "sliding_window_size" of Long, parameter
           "sliding_window_min_quality" of Long, parameter
           "leading_min_quality" of Long, parameter "trailing_min_quality" of
           Long, parameter "crop_length" of Long, parameter
           "head_crop_length" of Long, parameter "min_length" of Long
        :returns: instance of type "execTrimmomaticOutput" -> structure:
           parameter "output_filtered_ref" of type "data_obj_ref", parameter
           "output_unpaired_fwd_ref" of type "data_obj_ref", parameter
           "output_unpaired_rev_ref" of type "data_obj_ref", parameter
           "report" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN execTrimmomatic
        console = []
        self.log(console, 'Running execTrimmomatic with parameters: ')
        self.log(console, "\n"+pformat(input_params))
        report = ''
        trimmomatic_retVal = dict()
        trimmomatic_retVal['output_filtered_ref'] = None
        trimmomatic_retVal['output_unpaired_fwd_ref'] = None
        trimmomatic_retVal['output_unpaired_rev_ref'] = None

        token = ctx['token']
        wsClient = workspaceService(self.workspaceURL, token=token)
        headers = {'Authorization': 'OAuth '+token}
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token

        # object info
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        Set_types = ["KBaseSets.ReadsSet", "KBaseRNASeq.RNASeqSampleSet"]
        PE_types  = ["KBaseFile.PairedEndLibrary", "KBaseAssembly.PairedEndLibrary"]
        SE_types  = ["KBaseFile.SingleEndLibrary", "KBaseAssembly.SingleEndLibrary"]
        acceptable_types = Set_types + PE_types + SE_types

        # param checks
        required_params = ['input_reads_ref',
                           'output_ws',
                           'output_reads_name'
#                           'read_type'
                          ]
        for required_param in required_params:
            if required_param not in input_params or input_params[required_param] == None:
                raise ValueError ("Must define required param: '"+required_param+"'")

        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects']=[str(input_params['input_reads_ref'])]

        # Determine whether read library or read set is input object
        #
        try:
            input_reads_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':input_params['input_reads_ref']}]})[0]
            input_reads_obj_type = input_reads_obj_info[TYPE_I]
            input_reads_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_reads_obj_type)  # remove trailing version
            #input_reads_obj_version = input_reads_obj_info[VERSION_I]  # this is object version, not type version
        except Exception as e:
            raise ValueError('Unable to get read library object from workspace: (' + str(input_params['input_reads_ref']) +')' + str(e))

        if input_reads_obj_type not in acceptable_types:
            raise ValueError ("Input reads of type: '"+input_reads_obj_type+"'.  Must be one of "+", ".join(acceptable_types))

        # auto-detect reads type
        read_type = None
        if input_reads_obj_type in PE_types:
            read_type = 'PE'
        elif input_reads_obj_type in SE_types:
            read_type = 'SE'

        # get set
        #
        readsSet_ref_list = []
        readsSet_names_list = []
        if input_reads_obj_type in Set_types:
            try:
                #self.log (console, "INPUT_READS_REF: '"+input_params['input_reads_ref']+"'")  # DEBUG
                #setAPI_Client = SetAPI (url=self.callbackURL, token=ctx['token'])  # for SDK local.  doesn't work for SetAPI
                setAPI_Client = SetAPI (url=self.serviceWizardURL, token=ctx['token'], service_ver='beta')  # for dynamic service
                input_readsSet_obj = setAPI_Client.get_reads_set_v1 ({'ref':input_params['input_reads_ref'],'include_item_info':1})

            except Exception as e:
                raise ValueError('SetAPI FAILURE: Unable to get read library set object from workspace: (' + str(input_params['input_reads_ref'])+")\n" + str(e))
            for readsLibrary_obj in input_readsSet_obj['data']['items']:
                readsSet_ref_list.append(readsLibrary_obj['ref'])
                readsSet_names_list.append(readsLibrary_obj['info'][NAME_I])
                reads_item_type = readsLibrary_obj['info'][TYPE_I]
                reads_item_type = re.sub ('-[0-9]+\.[0-9]+$', "", reads_item_type)  # remove trailing version
                if reads_item_type in PE_types:
                    this_read_type = 'PE'
                elif reads_item_type in SE_types:
                    this_read_type = 'SE'
                else:
                    raise ValueError ("Can't handle read item type '"+reads_item_type+"' obj_name: '"+readsLibrary_obj['info'][NAME_I]+" in Set: '"+str(input_params['input_reads_ref'])+"'")
                if read_type != None and this_read_type != read_type:
                    raise ValueError ("Can't handle read Set: '"+str(input_params['input_reads_ref'])+"'.  Unable to process mixed PairedEndLibrary and SingleEndLibrary.  Please split into separate ReadSets")
                elif read_type == None:
                    read_type = this_read_type
        else:
            readsSet_ref_list = [input_params['input_reads_ref']]
            readsSet_names_list = [input_reads_obj_info[NAME_I]]


        # Iterate through readsLibrary members of set
        #
        report = ''
        trimmed_readsSet_ref = None
        unpaired_fwd_readsSet_ref = None
        unpaired_rev_readsSet_ref = None
        trimmed_RNASeqSampleSet_ref = None
        unpaired_fwd_SampleSet_ref = None
        unpaired_rev_SampleSet_ref = None

        trimmed_readsSet_refs      = []
        unpaired_fwd_readsSet_refs = []
        unpaired_rev_readsSet_refs = []

        for reads_item_i,input_reads_library_ref in enumerate(readsSet_ref_list):
            execTrimmomaticParams = { 'input_reads_ref': input_reads_library_ref,
                                      'output_ws': input_params['output_ws']
                                      }
            optional_params = [ #'read_type',
                               'adapterFa',
                               'seed_mismatches',
                               'palindrome_clip_threshold',
                               'simple_clip_threshold',
                               #'quality_encoding',
                               'translate_to_phred33',
                               'sliding_window_size',
                               'sliding_window_min_quality',
                               'leading_min_quality',
                               'trailing_min_quality',
                               'crop_length',
                               'head_crop_length',
                               'min_length'
                               ]
            for arg in optional_params:
                if arg in input_params:
                    execTrimmomaticParams[arg] = input_params[arg]

            # add auto-detected read_type
            execTrimmomaticParams['read_type'] = read_type

            # set output name
            if input_reads_obj_type not in Set_types:
                execTrimmomaticParams['output_reads_name'] = input_params['output_reads_name']
            else:
                execTrimmomaticParams['output_reads_name'] = readsSet_names_list[reads_item_i]+'_trimm'

            report += "RUNNING TRIMMOMATIC ON LIBRARY: "+str(input_reads_library_ref)+" "+str(readsSet_names_list[reads_item_i])+"\n"
            report += "-----------------------------------------------------------------------------------\n\n"

            # run Trimmomatic App for One Library at a Time
            trimmomaticSingleLibrary_retVal = self.execTrimmomaticSingleLibrary (ctx, execTrimmomaticParams)[0]

            # add to report
            report += trimmomaticSingleLibrary_retVal['report']+"\n\n"
            trimmed_readsSet_refs.append (trimmomaticSingleLibrary_retVal['output_filtered_ref'])
            unpaired_fwd_readsSet_refs.append (trimmomaticSingleLibrary_retVal['output_unpaired_fwd_ref'])
            unpaired_rev_readsSet_refs.append (trimmomaticSingleLibrary_retVal['output_unpaired_rev_ref'])


        # Just one Library
        if input_reads_obj_type not in ["KBaseSets.ReadsSet", "KBaseRNASeq.RNASeqSampleSet"]:

            # create return output object
            output = { 'report': report,
                       'output_filtered_ref': trimmed_readsSet_refs[0],
                       'output_unpaired_fwd_ref': unpaired_fwd_readsSet_refs[0],
                       'output_unpaired_rev_ref': unpaired_rev_readsSet_refs[0],
                     }
        # ReadsSet
        else:

            # save trimmed readsSet
            some_trimmed_output_created = False
            items = []
            for i,lib_ref in enumerate(trimmed_readsSet_refs):   # FIX: assumes order maintained
                if lib_ref == None:
                    #items.append(None)  # can't have 'None' items in ReadsSet
                    continue
                else:
                    some_trimmed_output_created = True
                    try:
                        label = input_readsSet_obj['data']['items'][i]['label']
                    except:
                        label = wsClient.get_object_info_new ({'objects':[{'ref':lib_ref}]})[0][NAME_I]
                    label = label + "_Trimm_paired"

                    items.append({'ref': lib_ref,
                                  'label': label
                                  #'data_attachment': ,
                                  #'info':
                                      })
            if some_trimmed_output_created:
                single_reads = False
                if read_type == 'SE':
                    reads_desc_ext = " Trimmomatic trimmed SingleEndLibrary"
                    reads_name_ext = "_trimm"
                    single_reads = True
                else:
                    reads_desc_ext = " Trimmomatic trimmed paired reads"
                    reads_name_ext = "_trimm_paired"
                output_readsSet_obj = { 'description': str(input_readsSet_obj['data']['description'])+reads_desc_ext,
                                        'items': items
                                        }
                output_readsSet_name = str(input_params['output_reads_name'])+reads_name_ext
                trimmed_readsSet_ref = setAPI_Client.save_reads_set_v1 ({'workspace_name': input_params['output_ws'],
                                                                         'output_object_name': output_readsSet_name,
                                                                         'data': output_readsSet_obj
                                                                         })['set_ref']
                trimmed_RNASeqSampleSet_ref = self._save_RNASeqSampleSet(
                                                            items,
                                                            input_params['output_ws'],
                                                            output_readsSet_name + '_SampleSet',
                                                            reads_desc_ext,
                                                            single_reads)
            else:
                self.log(console, "No trimmed output created")
                # raise ValueError ("No trimmed output created")


            # save unpaired forward readsSet
            some_unpaired_fwd_output_created = False
            if len(unpaired_fwd_readsSet_refs) > 0:
                items = []
                for i,lib_ref in enumerate(unpaired_fwd_readsSet_refs):  # FIX: assumes order maintained
                    if lib_ref == None:
                        #items.append(None)  # can't have 'None' items in ReadsSet
                        continue
                    else:
                        some_unpaired_fwd_output_created = True
                        try:
                            if len(unpaired_fwd_readsSet_refs) == len(input_readsSet_obj['data']['items']):
                                label = input_readsSet_obj['data']['items'][i]['label']
                            else:
                                label = wsClient.get_object_info_new ({'objects':[{'ref':lib_ref}]})[0][NAME_I]
                        except:
                            label = wsClient.get_object_info_new ({'objects':[{'ref':lib_ref}]})[0][NAME_I]
                        label = label + "_Trimm_unpaired_fwd"

                        items.append({'ref': lib_ref,
                                      'label': label
                                      #'data_attachment': ,
                                      #'info':
                                          })
                if some_unpaired_fwd_output_created:
                    output_readsSet_obj = { 'description': str(input_readsSet_obj['data']['description'])+" Trimmomatic unpaired fwd reads",
                                            'items': items
                                            }
                    output_readsSet_name = str(input_params['output_reads_name'])+'_trimm_unpaired_fwd'
                    unpaired_fwd_readsSet_ref = setAPI_Client.save_reads_set_v1 ({'workspace_name': input_params['output_ws'],
                                                                                  'output_object_name': output_readsSet_name,
                                                                                  'data': output_readsSet_obj
                                                                                  })['set_ref']
                    unpaired_fwd_SampleSet_ref = self._save_RNASeqSampleSet(
                                                            items,
                                                            input_params['output_ws'],
                                                            output_readsSet_name + '_SampleSet',
                                                            reads_desc_ext,
                                                            single_reads)
                else:
                    self.log (console, "no unpaired_fwd readsLibraries created")
                    unpaired_fwd_readsSet_ref = None

            # save unpaired reverse readsSet
            some_unpaired_rev_output_created = False
            if len(unpaired_rev_readsSet_refs) > 0:
                items = []
                for i,lib_ref in enumerate(unpaired_fwd_readsSet_refs):  # FIX: assumes order maintained
                    if lib_ref == None:
                        #item`s.append(None)  # can't have 'None' items in ReadsSet
                        continue
                    else:
                        some_unpaired_rev_output_created = True
                        try:
                            if len(unpaired_rev_readsSet_refs) == len(input_readsSet_obj['data']['items']):
                                label = input_readsSet_obj['data']['items'][i]['label']
                            else:
                                label = wsClient.get_object_info_new ({'objects':[{'ref':lib_ref}]})[0][NAME_I]

                        except:
                            label = wsClient.get_object_info_new ({'objects':[{'ref':lib_ref}]})[0][NAME_I]
                        label = label + "_Trimm_unpaired_rev"

                        items.append({'ref': lib_ref,
                                      'label': label
                                      #'data_attachment': ,
                                      #'info':
                                          })
                if some_unpaired_rev_output_created:
                    output_readsSet_obj = { 'description': str(input_readsSet_obj['data']['description'])+" Trimmomatic unpaired rev reads",
                                            'items': items
                                            }
                    output_readsSet_name = str(input_params['output_reads_name'])+'_trimm_unpaired_rev'
                    unpaired_rev_readsSet_ref = setAPI_Client.save_reads_set_v1 ({'workspace_name': input_params['output_ws'],
                                                                                  'output_object_name': output_readsSet_name,
                                                                                  'data': output_readsSet_obj
                                                                                  })['set_ref']

                    unpaired_rev_SampleSet_ref = self._save_RNASeqSampleSet(
                                                            items,
                                                            input_params['output_ws'],
                                                            output_readsSet_name + '_SampleSet',
                                                            reads_desc_ext,
                                                            single_reads)
                else:
                    self.log (console, "no unpaired_rev readsLibraries created")
                    unpaired_rev_readsSet_ref = None


            # create return output object
            output = {'report': report,
                      'output_filtered_ref': trimmed_readsSet_ref,
                      'output_unpaired_fwd_ref': unpaired_fwd_readsSet_ref,
                      'output_unpaired_rev_ref': unpaired_rev_readsSet_ref,
                      'output_filtered_sampleset_ref': trimmed_RNASeqSampleSet_ref,
                      'output_unpaired_sampleset_fwd_ref': unpaired_fwd_SampleSet_ref,
                      'output_unpaired_sampleset_rev_ref': unpaired_rev_SampleSet_ref
                    }

        #END execTrimmomatic

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method execTrimmomatic return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def execTrimmomaticSingleLibrary(self, ctx, input_params):
        """
        :param input_params: instance of type "execTrimmomaticInput"
           (execTrimmomatic() ** ** the local method that runs Trimmomatic on
           each read library) -> structure: parameter "input_reads_ref" of
           type "data_obj_ref", parameter "output_ws" of type
           "workspace_name" (** Common types), parameter "output_reads_name"
           of type "data_obj_name", parameter "read_type" of String,
           parameter "adapterFa" of String, parameter "seed_mismatches" of
           Long, parameter "palindrome_clip_threshold" of Long, parameter
           "simple_clip_threshold" of Long, parameter "translate_to_phred33"
           of type "bool", parameter "sliding_window_size" of Long, parameter
           "sliding_window_min_quality" of Long, parameter
           "leading_min_quality" of Long, parameter "trailing_min_quality" of
           Long, parameter "crop_length" of Long, parameter
           "head_crop_length" of Long, parameter "min_length" of Long
        :returns: instance of type "execTrimmomaticOutput" -> structure:
           parameter "output_filtered_ref" of type "data_obj_ref", parameter
           "output_unpaired_fwd_ref" of type "data_obj_ref", parameter
           "output_unpaired_rev_ref" of type "data_obj_ref", parameter
           "report" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN execTrimmomaticSingleLibrary
        console = []
        self.log(console, 'Running Trimmomatic with parameters: ')
        self.log(console, "\n"+pformat(input_params))
        report = ''
        retVal = dict()
        retVal['output_filtered_ref'] = None
        retVal['output_unpaired_fwd_ref'] = None
        retVal['output_unpaired_rev_ref'] = None

        token = ctx['token']
        wsClient = workspaceService(self.workspaceURL, token=token)
        headers = {'Authorization': 'OAuth '+token}
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token

        # object info
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        #Set_types = ["KBaseSets.ReadsSet", "KBaseRNASeq.RNASeqSampleSet"]
        PE_types = ["KBaseFile.PairedEndLibrary", "KBaseAssembly.PairedEndLibrary"]
        SE_types = ["KBaseFile.SingleEndLibrary", "KBaseAssembly.SingleEndLibrary"]
        acceptable_types = PE_types + SE_types

        # param checks
        required_params = ['input_reads_ref',
                           'output_ws',
                           'output_reads_name',
                           'read_type'
                          ]
        for required_param in required_params:
            if required_param not in input_params or input_params[required_param] == None:
                raise ValueError ("Must define required param: '"+required_param+"'")

        # and param defaults
        defaults = {
            #'quality_encoding':           'phred33',
            'seed_mismatches':            '0', # '2',
            'palindrome_clip_threshold':  '0', # '3',
            'simple_clip_threshold':      '0', # '10',
            'crop_length':                '0',
            'head_crop_length':           '0',
            'leading_min_quality':        '0', # '3',
            'trailing_min_quality':       '0', # '3',
            'sliding_window_size':        '0', # '4',
            'sliding_window_min_quality': '0', # '15',
            'min_length':                 '0', # '36'
        }
        for arg in defaults.keys():
            if arg not in input_params or input_params[arg] is None or input_params[arg] == '':
                input_params[arg] = defaults[arg]

        # conditional arg behavior
        arg = 'adapterFa'
        if arg not in input_params or input_params[arg] is None or input_params[arg] == '':
            input_params['adapterFa'] = None
            input_params['seed_mismatches'] = None
            input_params['palindrome_clip_threshold'] = None
            input_params['simple_clip_threshold'] = None


        #load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects']=[str(input_params['input_reads_ref'])]

        # Determine whether read library is of correct type
        #
        try:
            input_reads_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':input_params['input_reads_ref']}]})[0]
            input_reads_obj_type = input_reads_obj_info[TYPE_I]
            #input_reads_obj_version = input_reads_obj_info[VERSION_I]  # this is object version, not type version

        except Exception as e:
            raise ValueError('Unable to get read library object from workspace: (' + str(input_params['input_reads_ref']) +')' + str(e))

        input_reads_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_reads_obj_type)  # remove trailing version
        acceptable_types = PE_types + SE_types
        if input_reads_obj_type not in acceptable_types:
            raise ValueError ("Input reads of type: '"+input_reads_obj_type+"'.  Must be one of "+", ".join(acceptable_types))


        # Confirm user is paying attention (matters because Trimmomatic params are very different for PairedEndLibary and SingleEndLibrary
        #
        if input_params['read_type'] == 'PE' and not input_reads_obj_type in PE_types:
            raise ValueError ("read_type set to 'Paired End' but object is SingleEndLibrary")
        if input_params['read_type'] == 'SE' and not input_reads_obj_type in SE_types:
            raise ValueError ("read_type set to 'Single End' but object is PairedEndLibrary")


        # Instatiate ReadsUtils
        #
        try:
            readsUtils_Client = ReadsUtils (url=self.callbackURL, token=ctx['token'])  # SDK local

            readsLibrary = readsUtils_Client.download_reads ({'read_libraries': [input_params['input_reads_ref']],
                                                             'interleaved': 'false'
                                                             })
        except Exception as e:
            raise ValueError('Unable to get read library object from workspace: (' + str(input_params['input_reads_ref']) +")\n" + str(e))


        if input_params['read_type'] == 'PE':

            # Download reads Libs to FASTQ files
            input_fwd_file_path = readsLibrary['files'][input_params['input_reads_ref']]['files']['fwd']
            input_rev_file_path = readsLibrary['files'][input_params['input_reads_ref']]['files']['rev']
            sequencing_tech     = readsLibrary['files'][input_params['input_reads_ref']]['sequencing_tech']


            # DEBUG
#            self.log (console, "FWD_INPUT\n")
#            fwd_reads_handle = open (input_fwd_file_path, 'r')
#            for line_i in range(20):
#                self.log (console, fwd_reads_handle.readline())
#            fwd_reads_handle.close ()
#            self.log (console, "REV_INPUT\n")
#            rev_reads_handle = open (input_rev_file_path, 'r')
#            for line_i in range(20):
#                self.log (console, rev_reads_handle.readline())
#            rev_reads_handle.close ()


            # Set Params
            #
            trimmomatic_params  = self.parse_trimmomatic_steps(input_params)

            # add auto-detected quality_encoding
            if self.is_fastq_phred64 (input_fwd_file_path):
                quality_encoding = 'phred64'
            else:
                quality_encoding = 'phred33'

            trimmomatic_options = str(input_params['read_type']) + ' -' + quality_encoding
            self.log(console, pformat(trimmomatic_params))
            self.log(console, pformat(trimmomatic_options))


            # Run Trimmomatic
            #
            self.log(console, 'Starting Trimmomatic')
            input_fwd_file_path = re.sub ("\.fq$", "", input_fwd_file_path)
            input_fwd_file_path = re.sub ("\.FQ$", "", input_fwd_file_path)
            input_rev_file_path = re.sub ("\.fq$", "", input_rev_file_path)
            input_rev_file_path = re.sub ("\.FQ$", "", input_rev_file_path)
            input_fwd_file_path = re.sub ("\.fastq$", "", input_fwd_file_path)
            input_fwd_file_path = re.sub ("\.FASTQ$", "", input_fwd_file_path)
            input_rev_file_path = re.sub ("\.fastq$", "", input_rev_file_path)
            input_rev_file_path = re.sub ("\.FASTQ$", "", input_rev_file_path)
            output_fwd_paired_file_path   = input_fwd_file_path+"_trimm_fwd_paired.fastq"
            output_fwd_unpaired_file_path = input_fwd_file_path+"_trimm_fwd_unpaired.fastq"
            output_rev_paired_file_path   = input_rev_file_path+"_trimm_rev_paired.fastq"
            output_rev_unpaired_file_path = input_rev_file_path+"_trimm_rev_unpaired.fastq"
            input_fwd_file_path           = input_fwd_file_path+".fastq"
            input_rev_file_path           = input_rev_file_path+".fastq"

            cmdstring = " ".join( (self.TRIMMOMATIC, trimmomatic_options,
                            input_fwd_file_path,
                            input_rev_file_path,
                            output_fwd_paired_file_path,
                            output_fwd_unpaired_file_path,
                            output_rev_paired_file_path,
                            output_rev_unpaired_file_path,
                            trimmomatic_params) )

            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)

            outputlines = []
            while True:
                line = cmdProcess.stdout.readline()
                outputlines.append(line)
                if not line: break
                self.log(console, line.replace('\n', ''))

            cmdProcess.stdout.close()
            cmdProcess.wait()
            self.log(console, 'return code: ' + str(cmdProcess.returncode) + '\n')
            if cmdProcess.returncode != 0:
                raise ValueError('Error running kb_trimmomatic, return code: ' +
                                 str(cmdProcess.returncode) + '\n')

            #report += "cmdstring: " + cmdstring + " stdout: " + stdout + " stderr " + stderr

            # free up disk
            os.remove(input_fwd_file_path)
            os.remove(input_rev_file_path)

            # Only keep line that starts with Input and that
            # needs to be parsed for HTML report
            for line in outputlines:
                if line.startswith('Input'):
                    report += line

            # upload paired reads
            if not os.path.isfile (output_fwd_paired_file_path) \
                or os.path.getsize (output_fwd_paired_file_path) == 0 \
                or not os.path.isfile (output_rev_paired_file_path) \
                or os.path.getsize (output_rev_paired_file_path) == 0:
                retVal['output_filtered_ref'] = None
                report += "\n\nNo reads were trimmed, so no trimmed reads object was generated."
            else:

                # standardize quality encoding
                if 'translate_to_phred33' in input_params and input_params['translate_to_phred33'] == 1 and quality_encoding == 'phred64':
                #if False:  # DEBUG
                    self.log (console, "TRANSLATING OUTPUT FWD PAIRED FASTQ FILE...")
                    output_fwd_paired_file_path = self.translate_fastq_from_phred64_to_phred33 \
                                                  (output_fwd_paired_file_path, \
                                                   re.sub ("\.fastq$", ".q33.fastq", output_fwd_paired_file_path))
                    output_rev_paired_file_path = self.translate_fastq_from_phred64_to_phred33 \
                                                  (output_rev_paired_file_path, \
                                                   re.sub ("\.fastq$", ".q33.fastq", output_rev_paired_file_path))

                output_obj_name = input_params['output_reads_name']+'_paired'
                self.log(console, 'Uploading trimmed paired reads: '+output_obj_name)
                retVal['output_filtered_ref'] = readsUtils_Client.upload_reads ({ 'wsname': str(input_params['output_ws']),
                                                                                  'name': output_obj_name,
                                                                                  # remove sequencing_tech arg once ReadsUtils is updated to accept source_reads_ref
                                                                                  #'sequencing_tech': sequencing_tech,
                                                                                  'source_reads_ref': input_params['input_reads_ref'],
                                                                                  'fwd_file': output_fwd_paired_file_path,
                                                                                  'rev_file': output_rev_paired_file_path
                                                                                  })['obj_ref']

                # free up disk
                os.remove(output_fwd_paired_file_path)
                os.remove(output_rev_paired_file_path)


            # upload reads forward unpaired
            if not os.path.isfile (output_fwd_unpaired_file_path) \
                or os.path.getsize (output_fwd_unpaired_file_path) == 0:

                retVal['output_unpaired_fwd_ref'] = None
            else:
                # standardize quality encoding
                if 'translate_to_phred33' in input_params and input_params['translate_to_phred33'] == 1 and quality_encoding == 'phred64':
                #if False:  # DEBUG
                    self.log (console, "TRANSLATING OUTPUT FWD UNPAIRED FASTQ FILE...")
                    output_fwd_unpaired_file_path = self.translate_fastq_from_phred64_to_phred33 \
                                                    (output_fwd_unpaired_file_path, \
                                                     re.sub ("\.fastq$", ".q33.fastq", output_fwd_unpaired_file_path))

                output_obj_name = input_params['output_reads_name']+'_unpaired_fwd'
                self.log(console, '\nUploading trimmed unpaired forward reads: '+output_obj_name)
                retVal['output_unpaired_fwd_ref'] = readsUtils_Client.upload_reads ({ 'wsname': str(input_params['output_ws']),
                                                                                      'name': output_obj_name,
                                                                                      # remove sequencing_tech arg once ReadsUtils is updated to accept source_reads_ref
                                                                                      #'sequencing_tech': sequencing_tech,
                                                                                      'source_reads_ref': input_params['input_reads_ref'],
                                                                                      'fwd_file': output_fwd_unpaired_file_path
                                                                                      })['obj_ref']

                # free up disk
                os.remove(output_fwd_unpaired_file_path)

            # upload reads reverse unpaired
            if not os.path.isfile (output_rev_unpaired_file_path) \
                or os.path.getsize (output_rev_unpaired_file_path) == 0:

                retVal['output_unpaired_rev_ref'] = None
            else:
                # standardize quality encoding
                if 'translate_to_phred33' in input_params and input_params['translate_to_phred33'] == 1 and quality_encoding == 'phred64':
                #if False:  # DEBUG
                    self.log (console, "TRANSLATING OUTPUT REV UNPAIRED FASTQ FILE...")
                    output_rev_unpaired_file_path = self.translate_fastq_from_phred64_to_phred33 \
                                                    (output_rev_unpaired_file_path, \
                                                     re.sub ("\.fastq$", ".q33.fastq", output_rev_unpaired_file_path))

                output_obj_name = input_params['output_reads_name']+'_unpaired_rev'
                self.log(console, '\nUploading trimmed unpaired reverse reads: '+output_obj_name)
                retVal['output_unpaired_rev_ref'] = readsUtils_Client.upload_reads ({ 'wsname': str(input_params['output_ws']),
                                                                                      'name': output_obj_name,
                                                                                      # remove sequencing_tech arg once ReadsUtils is updated to accept source_reads_ref
                                                                                      #'sequencing_tech': sequencing_tech,
                                                                                      'source_reads_ref': input_params['input_reads_ref'],
                                                                                      'fwd_file': output_rev_unpaired_file_path
                                                                                      })['obj_ref']

                # free up disk
                os.remove(output_rev_unpaired_file_path)


        # SingleEndLibrary
        #
        else:
            self.log(console, "Downloading Single End reads file...")

            # Download reads Libs to FASTQ files
            input_fwd_file_path = readsLibrary['files'][input_params['input_reads_ref']]['files']['fwd']
            sequencing_tech     = readsLibrary['files'][input_params['input_reads_ref']]['sequencing_tech']

            # Set Params
            #
            trimmomatic_params  = self.parse_trimmomatic_steps(input_params)

            # add auto-detected quality_encoding
            if self.is_fastq_phred64 (input_fwd_file_path):
                quality_encoding = 'phred64'
            else:
                quality_encoding = 'phred33'

            trimmomatic_options = str(input_params['read_type']) + ' -' + quality_encoding
            self.log(console, pformat(trimmomatic_params))
            self.log(console, pformat(trimmomatic_options))


            # Run Trimmomatic
            #
            self.log(console, 'Starting Trimmomatic')
            input_fwd_file_path = re.sub ("\.fq$", "", input_fwd_file_path)
            input_fwd_file_path = re.sub ("\.FQ$", "", input_fwd_file_path)
            input_fwd_file_path = re.sub ("\.fastq$", "", input_fwd_file_path)
            input_fwd_file_path = re.sub ("\.FASTQ$", "", input_fwd_file_path)
            output_fwd_file_path = input_fwd_file_path+"_trimm_fwd.fastq"
            input_fwd_file_path  = input_fwd_file_path+".fastq"

            cmdstring = " ".join( (self.TRIMMOMATIC, trimmomatic_options,
                            input_fwd_file_path,
                            output_fwd_file_path,
                            trimmomatic_params) )

            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)

            #report += "cmdstring: " + cmdstring

            outputlines = []
            while True:
                line = cmdProcess.stdout.readline()
                outputlines.append(line)
                if not line: break
                self.log(console, line.replace('\n', ''))

            cmdProcess.stdout.close()
            cmdProcess.wait()
            self.log(console, 'return code: ' + str(cmdProcess.returncode) + '\n')
            if cmdProcess.returncode != 0:
                raise ValueError('Error running kb_trimmomatic, return code: ' +
                                 str(cmdProcess.returncode) + '\n')

            # Only keep line that starts with Input and that
            # needs to be parsed for HTML report

            for line in outputlines:
                if line.startswith('Input'):
                    report += line

            # free up disk
            os.remove(input_fwd_file_path)

            # get read count
            match = re.search(r'Surviving: (\d+)', report)
            readcount = match.group(1)

            # upload reads
            if not os.path.isfile (output_fwd_file_path) \
                or os.path.getsize (output_fwd_file_path) == 0:

                retVal['output_filtered_ref'] = None
            else:
                # standardize quality encoding
                if 'translate_to_phred33' in input_params and input_params['translate_to_phred33'] == 1 and quality_encoding == 'phred64':
                #if False:  # DEBUG
                    self.log (console, "TRANSLATING OUTPUT FASTQ FILE...")
                    output_fwd_file_path = self.translate_fastq_from_phred64_to_phred33 \
                                           (output_fwd_file_path, \
                                            re.sub ("\.fastq$", ".q33.fastq", output_fwd_file_path))

                output_obj_name = input_params['output_reads_name']
                self.log(console, 'Uploading trimmed reads: '+output_obj_name)

                retVal['output_filtered_ref'] = readsUtils_Client.upload_reads ({ 'wsname': str(input_params['output_ws']),
                                                                                  'name': output_obj_name,
                                                                                  # remove sequencing_tech arg once ReadsUtils is updated to accept source_reads_ref
                                                                                  #'sequencing_tech': sequencing_tech,
                                                                                  'source_reads_ref': input_params['input_reads_ref'],
                                                                                  'fwd_file': output_fwd_file_path
                                                                                  })['obj_ref']

                # free up disk
                os.remove(output_fwd_file_path)


        # return created objects
        #
        output = { 'report': report,
                   'output_filtered_ref': retVal['output_filtered_ref'],
                   'output_unpaired_fwd_ref': retVal['output_unpaired_fwd_ref'],
                   'output_unpaired_rev_ref': retVal['output_unpaired_rev_ref']
                 }
        #END execTrimmomaticSingleLibrary

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method execTrimmomaticSingleLibrary return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK", 'message': "", 'version': self.VERSION,
                     'git_url': self.GIT_URL, 'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
