import embed from 'vega-embed';

import buildCharts from './createChart.js';

buildCharts('data/Compare_MOTU.csv').then(chart => {
    embed("#vis", chart)
        // result.view provides access to the Vega View API
        .then(result => {
            // console.log(result.view.data("distances"));
        })
        .catch(console.warn);

});

