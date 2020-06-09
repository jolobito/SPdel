'use strict';
import * as d3Fetch from 'd3-fetch';

import dataParser from './dataParser.js';

function buildMatrixChart(){
    let matrixChart;
    let matrixMarks;
    let matrixScales;

    matrixScales = [
        {
            name: "position",
            type: "band",
            domain: {"data": "individuals", "field": "id", "sort": true},
            range: {"step": {"signal": "cellSize"}}
        },
        {
            name: "color",
            type: "ordinal",
            range: "category",
            domain: {
                "fields": [
                    {"data": "individuals", "field": "group"},
                ],
                sort: true
            }
        }
    ];

    matrixMarks = [
        {
            type: "rect",
            from: {"data": "distances"},
            encode: {
                enter: {
                    x: {"scale": "position", "field": "sourceNode.id"},
                    y: {"scale": "position", "field": "targetNode.id"},
                    width: {"scale": "position", "band": 1, "offset": -1},
                    height: {"scale": "position", "band": 1, "offset": -1},
                    fill: [
                        {test: "datum.sourceNode.group != datum.targetNode.group",
                            scale: "distance", "field": "distance"},
                        {scale: "distance_intra", "field": "distance"}
                        ],
                    tooltip:  {
                        signal: "{'there is a distance of ': datum.distance, " +
                            "'between': datum.sourceNode.id, " +
                            "'and': datum.targetNode.id}"
                    }
                },
                update: {
                    opacity:[
                        {test: " (!length(data('selected_dist'))  && !length(data('selected')) ) ||" +
                                "((  length(data('selected_dist')) &&" +
                                "  (data('selected_dist')[0].value.bin0 <= datum.distance  &&  " +
                                "   data('selected_dist')[0].value.bin1 > datum.distance ) " +
                                ") || " +
                                "(  length(data('selected')) && " +
                                " ( indata('selected', 'value', datum.sourceNode.group) && " +
                                "   indata('selected', 'value', datum.targetNode.group) ) " +
                                "))",
                            value: 0.9},
                        {"value": 0.15}
                    ],
                }
            }
        },
        {
            type: "rect",
            from: {"data": "individuals"},
            encode: {
                enter:{
                    tooltip:  {
                        signal: "{'title': datum.id, " +
                            "'group': datum.group," +
                            "'species': datum.species}"
                    }
                },
                update: {
                    "x": {"offset": -15},
                    "y": {"scale": "position", "field": "id", "offset": 2},
                    "width": {"scale": "position", "band": 0.5, "offset": -1},
                    "height": {"scale": "position", "band": 1, "offset": -1},
                    "fill": {"scale": "color", "field": "group"},
                    opacity:[
                        {"test": "!length(data('selected')) ||" +
                                "(indata('selected', 'value', datum.group) )",
                            value: 0.9},
                        {"value": 0.15}
                    ],
                }
            }
        },
        {
            type: "rect",
            from: {"data": "individuals"},
            encode: {
                enter:{
                    tooltip:  {
                        signal: "{'title': datum.id, " +
                            "'group': datum.group," +
                            "'species': datum.species}"
                    }
                },
                "update": {
                    "x": {"scale": "position", "field": "id", "offset": 2},
                    "y": {"offset": -15},
                    "width": {"scale": "position", "band": 1, "offset": -1},
                    "height": {"scale": "position", "band": 0.5, "offset": -1},
                    "fill": {"scale": "color", "field": "group"},
                    opacity:[
                        {"test": "!length(data('selected')) ||" +
                                "(indata('selected', 'value', datum.group) )",
                            value: 0.9},
                        {"value": 0.15}
                    ],
                }
            }
        },
        {
            type: "text",
            name: "columns",
            from: {"data": "individuals"},
            encode: {
                enter:{
                    tooltip:  {
                        signal: "{'title': datum.id, " +
                            "'group': datum.group," +
                            "'species': datum.species}"
                    }
                },
                "update": {
                    "x": {"scale": "position", "field": "id", "band": 0.5},
                    "y": {"offset": -25},
                    "text": {"field": "id"},
                    "fontSize": {"scale": "position", "band": 0.7},
                    "angle": {"value": -50},
                    "align": {"value": "left"},
                    "baseline": {"value": "middle"},
                    "fill": [
                        {"value": "black"}
                    ],
                    opacity:[
                        {"test": "!length(data('selected')) ||" +
                                "(indata('selected', 'value', datum.group) )",
                            value: 0.9},
                        {"value": 0.15}
                    ],
                }
            }
        },
        {
            type: "text",
            name: "rows",
            from: {"data": "individuals"},
            encode: {
                enter:{
                    tooltip:  {
                        signal: "{'title': datum.id, " +
                            "'group': datum.group," +
                            "'species': datum.species}"
                    }
                },
                "update": {
                    "x": {"offset": -25},
                    "y": {"scale": "position", "field": "id", "band": 0.8},
                    "text": {"field": "id"},
                    "fontSize": {"scale": "position", "band": 1},
                    "align": {"value": "right"},
                    "baseline": {"value": "middle"},
                    "fill": [
                        {"value": "black"}
                    ],
                    opacity:[
                        {"test": "!length(data('selected')) ||" +
                                "(indata('selected', 'value', datum.group) )",
                            value: 0.9},
                        {"value": 0.15}
                    ],
                }
            }
        },
    ];
    matrixChart = {
        description: "Distance matrix",
        name: "matrix",
        type: "group",
        encode: {
            enter: {
                y: {signal: "matrixY"},
                width: {signal: "width"},
                height: 600,
                fill: {value: "transparent"}
            }
        },
        scales: matrixScales,
        marks: matrixMarks
    };

    return matrixChart;
}

function histogramChart(name, description, data, pos, scale_color){
    let interHistogram;
    let scales;
    let axes;
    let marks;

    scales = [
        {
            name: "xscale",
            type: "linear",
            range: [0, {signal: "chartWidth"}],
            domain: {
                signal: "dist-extent"
            }
        },
        {
            name: "yscale",
            type: "linear",
            range: [{signal: "chartHeight"}, 0], round: true,
            domain: {"data": data, "field": "count"},
            zero: true, "nice": true
        }
    ];

    axes = [
        {"orient": "bottom", "scale": "xscale", "zindex": 1},
        {"orient": "left", "scale": "yscale", "tickCount": 5, "zindex": 1}
    ];

    marks = [
        {
            type: "rect",
            name: name+"Bar",
            from: {"data": data},
            encode: {
                enter:{
                    tooltip:  {
                        signal: "{'frequency': datum.count, " +
                            "'distance between': datum.bin0," +
                            "'and': datum.bin1}"
                    }
                },
                update: {
                    "x": {"scale": "xscale", "field": "bin0"},
                    "x2": {"scale": "xscale", "field": "bin1",
                        "offset": {"signal": "binStep > 0.02 ? -0.5 : 0"}},
                    "y": {"scale": "yscale", "field": "count"},
                    "y2": {"scale": "yscale", "value": 0},
                    "fill": {"scale": scale_color, "field": "bin0"},
                    "opacity": [
                        {"test": "(!length(data('selected_dist')) || " +
                                "indata('selected_dist', 'value.bin0', datum.bin0))",
                            "value": 0.7 },
                        {"value": 0.15}
                    ],
                },
                "hover": { "fill": {"value": "firebrick"} }
            }
        }
    ];

    interHistogram = {
        description: description,
        title: {
            text: description,
            anchor: "end",
            offset: -10,
            fontWeight: 200
        },
        name: name,
        type: "group",
        encode: {
            enter: {
                y:  {signal: pos},
                width: {signal: "width"},
                height: {value: 100},
                fill: {value: "transparent"}
            }
        },
        axes: axes,
        scales: scales,
        marks: marks
    };

    return interHistogram;
}

function minmaxChart(){

    let interHistogram;
    let scales;
    let axes;
    let marks;

    scales = [
        {
            name: "xscale",
            type: "linear",
            range: [0, {signal: "minMaxSize"}],
            domain: {
                signal: "dist-extent"
            }
        },
        {
            name: "yscale",
            type: "linear",
            range: [{signal: "minMaxSize"}, 0],
            domain: {
                signal: "dist-extent"
            }
        }
    ];

    axes = [
        {"orient": "bottom", "scale": "xscale", "zindex": 1},
        {"orient": "left", "scale": "yscale", "zindex": 1}
    ];

    marks = [
        {
            type: "symbol",
            name: "minmaxSymbol",
            // interactive: true,
            from: {"data": "min-max-distances"},
            encode: {
                enter: {
                    size: {value: 50},
                    tooltip:  {
                        signal: "{'group': datum['sourceNode\\.group'] }"
                    }
                },
                update: {
                    "x": {"scale": "xscale", "field": "max_intra"},
                    "y": {"scale": "yscale", "field": "min_inter"},
                    "fill": {"value": "steelblue"},
                    "opacity": [
                        {"test": "(!length(data('selected')) || " +
                                "indata('selected', 'value', datum['sourceNode\\.group']))",
                            "value": 0.7 },
                        {"value": 0.15}
                    ],
                    "zindex": {"value": 1}
                },
                hover: {
                    "fill": {"value": "firebrick"},
                    "fillOpacity": {"value": 1},
                    "zindex": {"value": 2}
                }
            }
        },
        {
            type: "rule",
            interactive: false,
            encode: {
                update: {
                    x: {value: 0},
                    x2: {signal: "minMaxSize"},
                    y: {signal: "minMaxSize"},
                    y2: {value: 0},
                    stroke: {value: "firebrick"},
                    strokeDash: {value: [8,8]},
                    zindex: {"value": 0}
                }
            }
        }
    ];

    interHistogram = {
        description: "Min Max Chart",
        title: {
            text: "min-max chart",
            anchor: "end",
            // offset: -10,
            fontWeight: 200
        },
        name: "minmax",
        type: "group",
        encode: {
            enter: {
                x:  {signal: "minmaxX"},
                y:  {signal: "minmaxY"},
                width: {signal: "minMaxSize"},
                height: {signal: "minMaxSize"},
                // fill: {value: "transparent"}
            }
        },
        axes: axes,
        scales: scales,
        marks: marks
    };

    return interHistogram;
}

async function buildCharts(fileName) {

    let marks;
    let dataSpec;
    let signals;
    let scales;

    const raw = await d3Fetch.text(fileName);
    const data = await dataParser(raw);

    dataSpec = [
        {
            name: "individuals",
            values: data.individuals,
            transform: [{
                type: "formula",
                as: "order",
                expr: "datum.group"
            }]
        },
        {
            name: "distances",
            values: data.distances,
            transform: [
                {
                    type: "lookup", "from": "individuals", "key": "id",
                    fields: ["i", "j"], "as": ["sourceNode", "targetNode"]
                },
                {
                    type: "formula", as: "group",
                    expr: "datum.sourceNode.group === datum.targetNode.group ? datum.sourceNode.group : count"
                },
                {
                    type: "extent",
                    field: "distance",
                    signal: "dist-extent"
                },
            ]
        },
        {
            name: "distances-intra",
            source: "distances",
            transform: [
                {
                    type: "filter",
                    expr: "datum.sourceNode.group == datum.targetNode.group"
                },
            ]
        },
        {
            name: "distances-inter",
            source: "distances",
            transform: [
                {
                    type: "filter",
                    expr: "datum.sourceNode.group != datum.targetNode.group"
                },
            ]
        },
        {
            name: "distances-inter-grouped",
            source: "distances-inter",
            transform: [
                {
                    type: "aggregate",
                    groupby: ["sourceNode.group"],
                    fields: ["distance"],
                    ops: ["min"], as: ["min_inter"]
                }
            ]
        },
        {
            name: "min-max-distances",
            source: "distances-intra",
            transform: [
                {
                    type: "aggregate",
                    groupby: ["sourceNode.group"],
                    fields: ["distance"],
                    ops: ["max"], as: ["max_intra"]
                },
                {
                    type: "lookup",
                    from: "distances-inter-grouped",
                    key: "sourceNode\\.group",
                    fields: ["sourceNode\\.group"],
                    values: ["min_inter"],
                    as: ["min_inter"]
                }
            ]
        },
        {
            name: "dist-binned-intra",
            source: "distances-intra",
            transform: [
                {
                    type: "bin", field: "distance",
                    extent: {signal: "dist-extent"},
                    step: {signal: "binStep"},
                },
                {
                    type: "aggregate",
                    key: "bin0", groupby: ["bin0", "bin1"],
                    fields: ["bin0"], ops: ["count"], as: ["count"]
                }
            ]
        },
        {
            name: "dist-binned-inter",
            source: "distances-inter",
            transform: [
                {
                    type: "bin", field: "distance",
                    extent: {signal: "dist-extent"},
                    step: {signal: "binStep"},
                },
                {
                    type: "aggregate",
                    key: "bin0", groupby: ["bin0", "bin1"],
                    fields: ["bin0"], ops: ["count"], as: ["count"]
                }
            ]
        },
        {
            "name": "selected",
            "on": [
                {"trigger": "clicked", "remove": true},
                {"trigger": "clicked", "insert": "clicked"},
                {"trigger": "clear", "remove": true},
            ]
        },
        {
            "name": "selected_dist",
            "on": [
                // {"trigger": "!shift", "remove": true},
                {"trigger": "interbar_clicked", "remove": true},
                {"trigger": "intrabar_clicked", "remove": true},
                {"trigger": "interbar_clicked", "insert": "interbar_clicked"},
                {"trigger": "intrabar_clicked", "insert": "intrabar_clicked"},
                {"trigger": "clear", "remove": true},
            ]
        }
    ];

    signals = [
        {name: "cellSize", value: 8},
        {name: "count", update: "length(data('individuals'))"},
        {name: "width", update: "span(range('position'))"},
        {name: "height", update: "width"},
        {name: "chartHeight", value: 100 },
        {name: "chartWidth", value: 500 },
        {name: "minMaxSize", value: 220 },
        {name: "matrixY", value: 400 },
        {name: "minmaxX", value: 600 },
        {name: "minmaxY", value: 0 },
        {name: "interHistogramY", value: 0 },
        {name: "intraHistogramY", value: 140 },
        {name: "binStep", value: 0.3},
        {
            "name": "clear", "value": true,
            "on": [
                {
                    "events": "view:mouseup[!event.item]",
                    "update": "true",
                    "force": true
                }
            ]
        },
        {
            "name": "clicked", "value": null,
            "on": [
                {
                    "events": "@minmaxSymbol:click",
                    "update": "{value: datum['sourceNode\\.group']}",
                    "force":  true
                }
            ]
        },
        {
            "name": "interbar_clicked", "value": null,
            "on": [
                {
                    "events": "@interBar:click",
                    "update": "{value: datum}",
                    "force":  true
                }
            ]
        },
        {
            "name": "intrabar_clicked", "value": null,
            "on": [
                {
                    "events": "@intraBar:click",
                    "update": "{value: datum}",
                    "force":  true
                }
            ]
        },
    ];

    scales = [
        {
            name: "distance",
            type: "linear",
            range: {scheme: "tealblues"},
            domain: {
                "fields": [
                    {"data": "distances", "field": "distance"},
                ],
            },
            // reverse: true
        },
        {
            name: "distance_intra",
            type: "linear",
            range: {scheme: "browns"},
            domain: {
                "fields": [
                    {"data": "distances", "field": "distance"},
                ],
            },
            // reverse: true
        }
     ];

    marks = [
        buildMatrixChart(),
        histogramChart(
            "intra",
            "histogram intra",
            "dist-binned-intra",
            "intraHistogramY",
            "distance_intra"),
        histogramChart("inter",
            "histogram inter",
            "dist-binned-inter",
            "interHistogramY",
            "distance"),
        minmaxChart()
    ];

    return {
        "$schema": "https://vega.github.io/schema/vega/v5.json",
        width: 1000,
        height: 1000,
        padding: 5,
        name: 'chart',
        signals: signals,
        data: dataSpec,
        scales: scales,
        marks: marks
    };
}

export default buildCharts;
