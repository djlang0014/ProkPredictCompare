Chart.register(ChartDataLabels);
//gets the info from CGI send over, converts to json for chart.js
var prodigalInfo = {{prodigalInfo|tojson|safe}};
var glimmerInfo = {{glimmerInfo|tojson|safe}};
var fgsInfo = {{fgsInfo|tojson|safe}};
var mgaInfo = {{mgaInfo|tojson|safe}};
var accession = "{{accession}}";

var data = {
        labels: ['% Genes Detected', '% Extra Predictions', '% Perfect Match', 'Median CDS Difference (NTs)', '% 3 Match', '% 5 Match', '% No Match'],
        datasets: [{
            label: 'Prodigal',
            data: prodigalInfo,
            backgroundColor: 'rgba(255, 99, 132, 0.2)', // Red
        }, {
            label: 'Glimmer',
            data: glimmerInfo,
            backgroundColor: 'rgba(75, 192, 192, 0.2)', // Green
        }, {
            label: 'FGS',
            data: fgsInfo,
            backgroundColor: 'rgba(255, 206, 86, 0.2)', // Yellow
        }, {
            label: 'MGA',
            data: mgaInfo,
            backgroundColor: 'rgba(153, 102, 255, 0.2)', // Purple
        }]
      };

      var ctx = document.getElementById('myChart').getContext('2d');
      var myChart = new Chart(ctx, {
          type: 'bar',
          data: data,
          options: {
              responsive: true,
              plugins: {
                  title: {
                      display: true,
                      text: accession + ' prediction metrics'
                  },
                  datalabels: {
                      color: '#000000',
                      display: true,  // always display labels
                      anchor: 'end',
                      align: 'top',
                      formatter: function(value, context) {
                          return value.toFixed(2);
                      }
                  }
              }
          }
      });
         