
import traceback
import re


from datarecord import Types

from cuberecord import OERecordComputeCube
from Standards import Fields

from orionclient.session import in_orion, OrionSession
from orionclient.types import File
from os import environ

from openeye import oedepict
import TrjAnalysisCubes.utils as utl

import base64


_clus_floe_report_header = """
<html>
<head>
<style>

  .cb-floe-report-wrapper * {
    box-sizing: border-box;
    font-family: 'Helvetica Neue', Helvetica;
  }

  .cb-floe-report__row {
    width:100%;
    display:flex;
  }
  .cb-floe-report__sidebar {
    width: 25%;
  }
  .cb-floe-report__content {
    width: 75%;
  }
  .cb-floe-report__column > * {
    width: 100%;
    height: auto;
  }
  .cb-floe-report__dev-alert {
    background-color: rgb(255, 80, 80);
    color: white;
    padding: .5em;
    font-weight: normal;
  }
  .cb-floe-report-element--analysis {
    /*white-space: pre-wrap;*/
    padding: 0 10px 20px 10px;
    font-size: calc(.5vh + .75vw + .25vmin);
  }

  div.cb-floe-report__analysis-table {
    width: 100%;
  }

  div.cb-floe-report__analysis-table-row {
    display: flex;
    justify-content: space-between;
    padding: 0 5px;
  }

  div.cb-floe-report__analysis-table-row:nth-child(1) {
    font-weight: bold;
    border-bottom: 1px solid black;
  }

  div.cb-floe-report__analysis-table-row > span:nth-child(1) {
    flex-basis: 20%;
    text-align: center;
  }

  div.cb-floe-report__analysis-table-row > span:nth-child(2) {
    flex-basis: 15%;
    text-align: right;
  }

  div.cb-floe-report__analysis-table-row > span:nth-child(3) {
    flex-basis: 50%;
    text-align: left;
    margin-left: auto;
  }

  h2.cb-floe-report-element--header {
    margin-bottom: 0;
    text-align: left;
  }
  h3.cb-floe-report-element--header {
    margin-bottom: 0;
    text-align: center;
  }
  .cb-floe-report__dev-alert small{
    color: white !important;
  }

  /* tab styles */
  div.cb-floe-report__tab-wrapper input {
    display: none;
  }

  div.cb-floe-report__tab-wrapper label {
    display: block;
    float: left;
    padding: 5px 10px;
    cursor: pointer;
    border: 2px solid black;
    margin-right: 2px;
  }

  div.cb-floe-report__tab-wrapper input:checked + label {
    background: black;
    color: white;
    cursor: default;
  }

  div.cb-floe-report__tab-wrapper div.cb-floe-report__tab-content {
    display: none;
    padding: 5px 10px;
    clear: left;
    border: 2px solid black;
  }
"""

_clus_floe_report_header2 = """
</style>

<script type="text/javascript">
  // fix2018aug24 document.addEventListener('DOMContentLoaded', function() {
    var svgs = document.querySelectorAll('.cb-floe-report__tab-content svg');

    // Trigger a click on the buttons in all other SVGs that have a corresponding index
    function syncOtherButtons(svgs, idxSVG, idxBtn) {
      svgs.forEach(function(element, idxA, listA) {
        if (idxA == idxSVG) return false;
        var buttons = element.querySelectorAll('g[id^="tab_button_"]');
        buttons.forEach(function(button, idxB, listB) {
          if (idxB !== idxBtn) return false;
          var evt = new Event('click');
          evt.detail = evt.detail || {};
          evt.detail.isSynthetic = true;
          button.dispatchEvent(evt);
        });
      });
    };

    function listenForClick(e) {
      e.detail = e.detail || {};
      // use a flag to prevent scripted click events from causing infinite recursion!
      if (e.detail && !e.detail.isSynthetic) {
        syncOtherButtons(svgs, idxSVG, idxBtn);
        e.stopPropagation();
      }
    }

    // When a button inside an SVG is clicked, trigger clicks on the others
    svgs.forEach(function(element, idxSVG, listA) {
      var buttons = element.querySelectorAll('g[id^="tab_button_"]');
      buttons.forEach(function(button, idxBtn, listB) {
        button.addEventListener('click', function(e) {
          e.detail = e.detail || {};
          // use a flag to prevent scripted click events from causing infinite recursion!
          if (e.detail && !e.detail.isSynthetic) {
            syncOtherButtons(svgs, idxSVG, idxBtn);
            e.stopPropagation();
          }
        });
      });
    });
  // fix2018aug24 });
</script>
"""

_clus_floe_report_midHtml0 = """

</head>
<body>
<div class="cb-floe-report-wrapper">
  <h3 class="cb-floe-report__dev-alert">
    NOTE: this is an Alpha Test version.
    <br>
    <small>We are actively working on improving the MD sampling.</small>
  </h3>
  <div class="cb-floe-report__row">
  <div class="cb-floe-report__column cb-floe-report__sidebar">
    {query_depiction}
"""

_clus_floe_report_midHtml1 = """
  </div>

  <div class="cb-floe-report__column cb-floe-report__content">
    <h2> Analysis of Short Trajectory MD </h2>

    <div class="cb-floe-report__tab-wrapper">

"""

_clus_floe_report_midHtml2 = """      </div>
    </div>
  </div>

  <div class="cb-floe-report__row">
    <h2 style="text-align: center; width: 100%"> Ligand Clustering based on Active Site Alignment </h2>
  </div>

  <div class="cb-floe-report__row">
    <div class="cb-floe-report__column cb-floe-report__sidebar">
"""

_clus_floe_report_Trailer = """    </div>

    <div class="cb-floe-report__column cb-floe-report__content">
      <h3 class="cb-floe-report-element--header"> Cluster membership of ligand by Trajectory frame </h3>
      {clusters}
      <h3 class="cb-floe-report-element--header"> RMSD of ligand compared to initial pose, colored by cluster </h3>
      {rmsdInit}
    </div>
  </div>

</div>

</body>
</html>"""


def MakeClusterInfoText(dataDict, rgbVec):
    # Generate text string about Clustering information
    #
    text = []
    nFrames = dataDict['nFrames']
    text.append("""
      <br/><br/>
      <div class="cb-floe-report-element--analysis">""")
    text.append('Clustering by ligand RMSD after alignment by active site C_alphas:\n' )
    text.append('        <br>- Cluster method {}\n'.format( dataDict['ClusterMethod']) )
    text.append('        <br>- Using alpha={:.2f}\n'.format( dataDict['HDBSCAN_alpha']))
    text.append('        <br>- Clustered {} frames\n'.format(nFrames) )
    #
    if dataDict['nClusters']<2:
        text.append('        <br>- Produced {} cluster'.format( dataDict['nClusters']))
    else:
        text.append('        <br>- Produced {} clusters'.format( dataDict['nClusters']))
    nOutliers = dataDict['ClusterVec'].count(-1)
    text.append(' with {:4d} outliers\n'.format( nOutliers))
    text.append('<br><br>Clusters larger than 2%: ' )
    #
    text.append("""
        <div class="cb-floe-report__analysis-table">
          <div class="cb-floe-report__analysis-table-row">
            <span>Cluster</span>
            <span>Size</span>
            <span>Status</span>
          </div>\n
""")
    #
    for i, (count,rgb) in enumerate(zip(dataDict['ClusterCounts'],rgbVec)):
        percent = int(100*count/nFrames)
        if percent<2:
            continue
        status = 'major'
        if percent<10:
            status = 'minor (no analysis)'
        text.append("""
          <div class="cb-floe-report__analysis-table-row" style="
                background-color: rgb({r}, {g}, {b} );
                color: white;">
            <span>{clusID}</span>
            <span>{percent}%</span>
            <span>{status}</span>
          </div>\n""".format( clusID=i, percent=percent, status=status,
                            r=rgb[0], g=rgb[1], b=rgb[2]))
    #
    text.append("""
        </div>
      </div>
    """)
    #
    return text


def png_to_data_url(png_data):
    return "<img src='data:image/png;base64," + base64.b64encode(png_data).decode('utf-8') + "'>"


def trim_svg(svg):
    run_init = "\nif (init != null) { try {init() } catch(error) {  true; } };\n"
    svg = re.sub('(<script[^>]*> *<!\[CDATA\[)', '\\1' + run_init, svg)

    idx = svg.find('<svg')
    return svg[idx:]


#class MDTrajAnalysisClusterReport(ParallelMixin, OERecordComputeCube):
class MDTrajAnalysisClusterReport(OERecordComputeCube):
    title = 'Extract relevant outputs of MD Traj Cluster  Analysis'

    version = "0.1.0"
    classification = [["Simulation", "Traj Analysis"]]
    tags = ['Parallel Cube']

    description = """
    Extract relevant outputs of Ligand and Protein
    Short Traj MD Traj Analysis and write them to files.

    This cube takes as input the OERecord containing the work
    product of trajectory analysis on Short Traj MD results.
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        return

    def process(self, record, port):
        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)

            # title of entire solvated protein-ligand system
            opt['Logger'].info('Starting Floe Report generation for MD Traj Analysis')
            system_title = utl.RequestOEFieldType( record, Fields.title)
            opt['Logger'].info('{} Attempting to extract MD Traj Analysis results'
                .format(system_title) )
            ligInitPose = utl.RequestOEFieldType( record, Fields.ligand)
            protInitPose = utl.RequestOEFieldType( record, Fields.protein)
            asiteSVG = utl.PoseInteractionsSVG( ligInitPose, protInitPose, width=400, height=265)

            # Extract the traj SVG from the OETraj record
            analysesDone = utl.RequestOEField( record, 'AnalysesDone', Types.StringVec)
            if 'OETraj' not in analysesDone:
                raise ValueError('{} does not have OETraj analyses done'.format(system_title) )
            else:
                opt['Logger'].info('{} found OETraj analyses'.format(system_title) )
            # Extract the relevant traj SVG from the OETraj record
            oetrajRecord = utl.RequestOEField( record, 'OETraj', Types.Record)
            opt['Logger'].info('{} found OETraj record'.format(system_title) )
            trajSVG = utl.RequestOEField( oetrajRecord, 'TrajSVG', Types.String)

            # Extract the three plots from the TrajClus record
            analysesDone = utl.RequestOEField( record, 'AnalysesDone', Types.StringVec)
            if 'TrajClus' not in analysesDone:
                raise ValueError('{} does not have TrajClus analyses done'.format(system_title) )
            else:
                opt['Logger'].info('{} found TrajClus analyses'.format(system_title) )
            # Extract the relevant traj SVG from the TrajClus record
            clusRecord = utl.RequestOEField( record, 'TrajClus', Types.Record)
            opt['Logger'].info('{} found TrajClus record'.format(system_title) )
            trajHistRMSD_svg = utl.RequestOEField( clusRecord, 'HistSVG', Types.String)
            trajClus_svg = utl.RequestOEField( clusRecord, 'ClusSVG', Types.String)
            rmsdInit_svg = utl.RequestOEField( clusRecord, 'rmsdInitPose', Types.String)
            clusTrajSVG = utl.RequestOEField( clusRecord, 'ClusTrajSVG', Types.StringVec)
            opt['Logger'].info('{} found the TrajClus plots'.format(system_title) )

            # Generate text string about Clustering information
            clusData = {}
            clusData['nFrames'] = utl.RequestOEField(clusRecord, 'nFrames', Types.Int)
            clusData['ClusterMethod'] = utl.RequestOEField(clusRecord, 'ClusterMethod', Types.String)
            clusData['HDBSCAN_alpha'] = utl.RequestOEField(clusRecord, 'HDBSCAN_alpha', Types.Float)
            clusData['nClusters'] = utl.RequestOEField(clusRecord, 'nClusters', Types.Int)
            clusData['ClusterVec'] = utl.RequestOEField( clusRecord, 'Clusters', Types.IntVec)
            clusData['ClusterCounts'] = utl.RequestOEField( clusRecord, 'ClusterCounts', Types.IntVec)

            opt['Logger'].info('{} finished writing analysis files'.format(system_title) )

            # prepare the 2D structure depiction
            oedepict.OEPrepareDepiction(ligInitPose)
            img = oedepict.OEImage(400, 300)
            oedepict.OERenderMolecule(img, ligInitPose)

            # get the palette of graph marker colors
            nClustersP1 = clusData['nClusters']+1
            clusRGB = utl.ColorblindRGBMarkerColors( nClustersP1)
            clusRGB[-1] = (76, 76, 76)

            # write the report
            reportFName = system_title+'_ClusReport.html'
            report_file = open( reportFName, 'w')

            report_file.write(_clus_floe_report_header)

            for i in range(len(clusTrajSVG)+2):
                report_file.write("""
              div.cb-floe-report__tab-wrapper input:nth-of-type({clusID}):checked ~ .cb-floe-report__tab-content:nth-of-type({clusID}) {{ display: block; }}
            """.format( clusID=i+1))

            report_file.write(_clus_floe_report_header2)

            report_file.write(_clus_floe_report_midHtml0.format(
                query_depiction=oedepict.OEWriteImageToString("svg", img).decode("utf8")))

            analysis_txt = MakeClusterInfoText( clusData,clusRGB)
            report_file.write("".join(analysis_txt))

            report_file.write(_clus_floe_report_midHtml1 )

            report_file.write("""      <input type="radio" name="tab" id="cb-floe-report__tab-1-header" checked>
                  <label class="cb-floe-report__tab-label" for="cb-floe-report__tab-1-header">Overall</label>

            """)
            CurrentTabId = 1
            for i, (clus,rgb) in enumerate(zip(clusTrajSVG,clusRGB)):
                CurrentTabId = i+2
                report_file.write("""      <input type="radio" name="tab" id="cb-floe-report__tab-{tabID}-header">
                  <label class="cb-floe-report__tab-label" for="cb-floe-report__tab-{tabID}-header" style="
                            background-color: rgb({r},{g},{b});
                            color: white;">Cluster {clusNum}</label>

            """.format( tabID=CurrentTabId, clusNum=i, r=rgb[0], g=rgb[1], b=rgb[2]))
            report_file.write("""      <input type="radio" name="tab" id="cb-floe-report__tab-{tabID}-header">
                  <label class="cb-floe-report__tab-label" for="cb-floe-report__tab-{tabID}-header">Initial Pose</label>

            """.format( tabID=CurrentTabId+1, clusNum=i, r=rgb[0], g=rgb[1], b=rgb[2]))

            report_file.write("""      <div class="cb-floe-report__tab-content">
                    {traj}
                  </div>
            """.format(traj=trim_svg(trajSVG)) )
            for clusSVG in clusTrajSVG:
                report_file.write("""      <div class="cb-floe-report__tab-content">
                    {traj}
                  </div>
            """.format(traj=trim_svg(clusSVG)) )
            report_file.write("""      <div class="cb-floe-report__tab-content">
                    {traj}
                  </div>
            """.format(traj=trim_svg(asiteSVG)) )

            report_file.write(_clus_floe_report_midHtml2)

            report_file.write(_clus_floe_report_Trailer.format(
                clusters=trim_svg(trajClus_svg),
                rmsdInit=trim_svg(rmsdInit_svg)))

            report_file.close()

            if in_orion():
                session = OrionSession()

                ligName = ligInitPose.GetTitle()
                id_plmd = utl.RequestOEFieldType( record, Fields.id)
                file_upload = File.upload(session, "Report{} {}".format(id_plmd, ligName), "./"+reportFName)
                session.tag_resource(file_upload, "floe_report")
                job_id = environ.get('ORION_JOB_ID')
                if job_id:
                    session.tag_resource(file_upload, "Job {}".format(job_id))

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return

