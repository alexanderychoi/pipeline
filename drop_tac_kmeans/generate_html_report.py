#######################################
## generate_html_report.py allows you 
## generate an html report for alignment
## statistics 
#######################################
## Developed by:
## Paul Rivaud & Graham Heimberg
## paulrivaud.info@gmail.com
## gheimberg@gmail.com
## 2015
#######################################

def generate_html_report(common_path, name, bowtie_score, reads_counted, total_reads, preprocessing_saved_reads, 
						 dismissed_reads, dis_tso, dis_no_tac, dis_redund, dict_quality, dict_quality_scores):
	html_file = open(common_path+name+'_report.html', 'w+')
	html_file.write('''<!DOCTYPE html>
	<html>
		<head>
			<link rel="stylesheet" type="text/css" href="https://bootswatch.com/cerulean/bootstrap.min.css">
		</head>
		<body>
			<nav class="navbar navbar-default">
			  <div class="container-fluid">
			    <div class="navbar-header">
			      <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#bs-example-navbar-collapse-1">
			        <span class="sr-only">Toggle navigation</span>
			        <span class="icon-bar"></span>
			        <span class="icon-bar"></span>
			        <span class="icon-bar"></span>
			      </button>
			      <a class="navbar-brand" onclick="return false">HTML Report for ''')
	html_file.write(name)
	html_file.write('''</a>
			    </div>
			    <div class="collapse navbar-collapse" id="bs-example-navbar-collapse-1">
			      	<ul class="nav navbar-nav navbar-right">
		        		<li><a onclick="return false">Thomson Lab</a></li>
		      		</ul>
			    </div>
			  </div>
			</nav>
			<br>
			<h1 class="text-primary">Preprocessing results</h1>
			<br>
			<div id="preprocessing_bar" style="height: 300px; width: 50%;"></div>
			<br>
			<h1 class="text-primary">Dismissed reads information</h1>
			<br>
			<div id="dismissed_info" style="height: 300px; width: 50%;"></div>
			<br>
			<h1 class="text-primary">Bowtie2 overall alignment rate</h1>
			<br>
			<h1>''')
	html_file.write(str(bowtie_score))
	html_file.write('''%</h1>
			<br>
			<h1 class="text-primary">Reads used for gene expression quantification</h1>
			<br>
			<h1>''')
	html_file.write(str(round(100.0*reads_counted/total_reads,1)))
	html_file.write('''%</h1>
			<br>
			<h1 class="text-primary">Reads per cell distribution</h1>
			<br>
			<div id="readsDistribution" style="height: 300px; width: 50%;"></div>
			<br>
			<h1 class="text-primary">Quality scores histogram</h1>
			<br>
			<div id="qualityScore" style="height: 300px; width: 50%;"></div>
			


			</body>
		</html>


		<script type="text/javascript" src="http://canvasjs.com/assets/script/canvasjs.min.js"></script>
		<script type="text/javascript">
			function preprocessingPlot() {
				var chart = new CanvasJS.Chart("preprocessing_bar", {
							theme: "theme2",//theme1
							title:{
								text: ""
							},
							animationEnabled: true,   // change to true
							data: [              
							{
								// Change type to "bar", "splineArea", "area", "spline", "pie",etc.
								type: "column",
								dataPoints: [
									{ label: "Total reads",  y: ''')
	html_file.write(str(total_reads))
	html_file.write('''  },
						{ label: "Saved reads", y: ''')
	html_file.write(str(preprocessing_saved_reads))
	html_file.write('''  },
						{ label: "Dismissed reads", y: ''')
	html_file.write(str(dismissed_reads))
	html_file.write('''  }
								]
							}
							]
						});
						chart.render();
			}

			function dismissedInfo(){
				var chart1 = new CanvasJS.Chart("dismissed_info",
							{
								title:{
									text: ""
								},
					                        animationEnabled: true,
								theme: "theme2",
								data: [
								{        
									type: "doughnut",
									indexLabelFontFamily: "Garamond",       
									indexLabelFontSize: 20,
									startAngle:0,
									indexLabelFontColor: "dimgrey",       
									indexLabelLineColor: "darkgrey", 
					

									dataPoints: [
									{  y: ''')
	html_file.write(str(dis_tso))
	html_file.write(''', label: "Contains TSO" },
						{  y: ''')
	html_file.write(str(dis_no_tac))
	html_file.write(''', label: "No TAC" },
						{  y: ''')
	html_file.write(str(dis_redund))
	html_file.write(''', label: "Redundant" }

									]
								}
								]
							});
							chart1.render();
			}

			function distribution(){
				var chart2 = new CanvasJS.Chart("readsDistribution",
				    {
				      title:{
				      text: ""   
				      },
				      axisY:{
				        title:"Number of reads"   
				      },
				      animationEnabled: true,
				      data: [
				      {        
				        type: "stackedColumn",
				        toolTipContent: "{label}<br/><span style='\\"'color: {color};'\\"'><strong>{name}</strong></span>: {y} reads",
				        name: "Low quality",
				        showInLegend: "true",
				        dataPoints: [''')
	for key in sorted(dict_quality.keys()):
		html_file.write('{  y: ')
		html_file.write(str(dict_quality[key]['low']))
		html_file.write(', label:"')
		html_file.write(key)
		html_file.write('"},\n')
	html_file.write(''']

				      },  {        
				        type: "stackedColumn",
				        toolTipContent: "{label}<br/><span style='\\"'color: {color};'\\"'><strong>{name}</strong></span>: {y} reads",
				        name: "Good quality",
				        showInLegend: "true",
				        dataPoints: [''')
	for key in sorted(dict_quality.keys()):
		html_file.write('{  y: ')
		html_file.write(str(dict_quality[key]['high']))
		html_file.write(', label:"')
		html_file.write(key)
		html_file.write('"},\n')
	html_file.write(''']
				      }            
				      ]
				      ,
				      legend:{
				        cursor:"pointer",
				        itemclick: function(e) {
				          if (typeof (e.dataSeries.visible) ===  "undefined" || e.dataSeries.visible) {
					          e.dataSeries.visible = false;
				          }
				          else
				          {
				            e.dataSeries.visible = true;
				          }
				          chart2.render();
				        }
				      }
				    });

				    chart2.render();
			}

			function quality(){
				var chart3 = new CanvasJS.Chart("qualityScore",
					{
						animationEnabled: true,
						title:{
							text: ""
						},
						data: [
						{
							type: "column", //change type to bar, line, area, pie, etc
							dataPoints: [''')
	for key in dict_quality_scores.keys():
		html_file.write('{  x: ')
		html_file.write(key)
		html_file.write(', y: Math.log10(')
		html_file.write(str(dict_quality_scores[key]))
		html_file.write(')},\n')
	html_file.write(''']
						}
						]
					});

					chart3.render();
			}
			function loadAll() {
				preprocessingPlot();
				dismissedInfo();
				distribution();
				quality();
			}

			window.onload = loadAll;
		</script>''')
	html_file.close()
	print "HTML report completed"	
