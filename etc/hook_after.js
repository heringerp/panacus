/* global bootstrap: false */
(() => {
'use strict'
const tooltipTriggerList = Array.from(document.querySelectorAll('[data-bs-toggle="tooltip"]'))
tooltipTriggerList.forEach(tooltipTriggerEl => {
new bootstrap.Tooltip(tooltipTriggerEl)
})
})()

const pluginCanvasBackgroundColor = {
  id: 'customCanvasBackgroundColor',
  beforeDraw: (chart, args, options) => {
    const {ctx, chartArea: { top, bottom, left, right, width, height },
        scales: {x, y}
    } = chart;
    ctx.save();
    ctx.globalCompositeOperation = 'destination-over';
    ctx.fillStyle = options.color || '#99ffff';
    ctx.fillRect(left, top, width, height);
    ctx.restore();
  }
}

// Adapted from https://github.com/vega/vega-embed/
function post_to_vega_editor(window, data) {
    const url = 'https://vega.github.io/editor/';
    const editor = window.open(url);
    const wait = 10000;
    const step = 250;
    const {origin} = new URL(url);

    let count = ~~(wait / step);

    function listen(evt) {
        if (evt.source === editor) {
            count = 0;
            window.removeEventListener('message', listen, false);
        }
    }
    window.addEventListener('message', listen, false);

    function send() {
        if (count <= 0) {
            return;
        }
        editor.postMessage(data, origin);
        setTimeout(send, step);
        count -= 1;
    }
    setTimeout(send, step);
}

document.getElementById('btn-download-config').onclick = function() {
    let blob = new Blob([objects.config.first], {type: 'text/plain'});
    var a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = 'config.yaml';
    a.click();
}

for (let key in objects.datasets) {
    let element = objects.datasets[key];
    if (element instanceof Bar) {
        let h = element;
        let ctx = document.getElementById('chart-bar-' + h.id);
        let id = 'chart-bar-' + h.id;
        let data = {};
        if (h.ordinal) {
            data.values = h.data.values.map(d => ({
                ...d,
                label: Number(d.label)
            }));
        } else {
            data.values = h.data.values;
        }
        if (h.data.values.length >= 200 && h.ordinal) {
            mark_type = "area";
            x_encoding = {
                "field": "label",
                "title": h.x_label,
                "type": "quantitative",
                "scale": {
                    "nice": false
                }
            };
        } else {
            mark_type = "bar";
            x_encoding = {
                "field": "label",
                "title": h.x_label,
                "sort": null,
                "type": h.ordinal ? 'ordinal' : 'nominal',
                "axis": {
                    "labelOverlap": "greedy"
                }
            }
        }
        if (h.log_toggle) {
            data.values = data.values.filter((el) => el.value > 0);
        }
        let yourVlSpec = {
            $schema: 'https://vega.github.io/schema/vega-lite/v6.json',
            description: 'Bar',
            width: 1000,
            "autosize": {
                "type": "fit",
                "contains": "padding"
            },
            height: 350,
            // data: h.data,
            data,
            layer: [
                {
                    "mark": {"type": mark_type, "tooltip": true},
                    "encoding": {
                        "x": x_encoding,
                        "y": {"field": 'value', "title": h.y_label},
                    },
                },
            ]
        };

        function render(scaleType, thisId, vlSpec, add_listeners) {
            const copied_spec = JSON.parse(JSON.stringify(vlSpec)); // deep copy
            if (scaleType == "log") {
                copied_spec.layer[0].encoding.y.scale = { type: "log", domainMin: 1 }; // set scale type
                copied_spec.layer[0].encoding.y2 = { datum: 1 }; // set scale type
            } else {
                copied_spec.layer[0].encoding.y.scale = { type: "linear" }; // set scale type
                if ('y2' in copied_spec.layer[0].encoding) {
                    delete copied_spec.layer[0].encoding[y2]; // set scale type
                }
            }
            let opt = {
                "actions": false,
            };
            vegaEmbed(`#${CSS.escape(thisId)}`, copied_spec, opt).then(({ view, spec, vgSpec }) => {
                if (add_listeners) {
                    // Export PNG
                    let png_button = document.getElementById('btn-download-plot-png-' + h.id);
                    png_button.addEventListener('click', () => {
                        view.toImageURL('png').then(url => {
                            const a = document.createElement('a');
                            a.href = url;
                            a.download = 'visualization.png';
                            a.click();
                        });
                    });

                    // Export SVG
                    let svg_button = document.getElementById('btn-download-plot-svg-' + h.id);
                    svg_button.removeEventListener('click', svg_button);
                    svg_button.addEventListener('click', function svg_button() {
                        view.toImageURL('svg').then(url => {
                            const a = document.createElement('a');
                            a.href = url;
                            a.download = 'visualization.svg';
                            a.click();
                        });
                    });

                    // Open in Vega Editor
                    let vega_editor_button = document.getElementById('btn-download-plot-vega-editor-' + h.id);
                    vega_editor_button.addEventListener('click', () => {
                        post_to_vega_editor(window, {
                            mode: 'vega-lite',
                            spec: JSON.stringify(spec, null, 2),
                            renderer: undefined,
                            config: undefined,
                        });
                    });
                }
            });
        }

        if (h.log_toggle) {
            document.getElementById('btn-logscale-plot-bar-' + h.id).addEventListener('change', (event) => {
                if (event.currentTarget.checked) {
                    render("log", id, yourVlSpec, false);
                } else {
                    render("linear", id, yourVlSpec, false);
                }
            });
        }

        if (document.getElementById('btn-logscale-plot-bar-' + h.id).checked) {
            render("log", id, yourVlSpec, true);
        } else {
            render("linear", id, yourVlSpec, true);
        }
    } else if (element instanceof MultiBar) {
        let m = element;
        var ctx = document.getElementById('chart-bar-' + m.id);
        let id = 'chart-bar-' + m.id;
        let mark_type, x_encoding;
        if (m.data.values.length >= 300 && m.ordinal) {
            mark_type = "area";
            x_encoding = {
                "field": "label",
                "title": m.x_label,
                "type": "quantitative",
                "scale": {
                    "nice": false
                }
            };
        } else if (m.data.values.length >= 300) {
            mark_type = "area";
            x_encoding = {
                "field": "label",
                "title": m.x_label,
                "type": 'nominal',
                "sort": null,
                "axis": {
                    "labelOverlap": "greedy"
                }
            }
        } else {
            mark_type = "bar";
            x_encoding = {
                "field": "label",
                "title": m.x_label,
                "type": m.ordinal ? 'ordinal' : 'nominal',
                "sort": null,
                "axis": {
                    "labelOverlap": "greedy"
                }
            }
        }
        let layer_values = [
                {
                    "params": [
                        {
                          "name": "toggle-curves",
                          "bind": "legend",
                          "select": {
                            "type": "point",
                            "toggle": "true",
                            "fields": ["name"]
                          }
                        }
                    ],
                    "mark": {"type": mark_type, "tooltip": {"content": "data"}},
                    "encoding": {
                        "x": x_encoding,
                        "y": {
                            "field": "value",
                            "title": m.y_label,
                            "stack": null
                        },
                        "color": {
                            "field": "name",
                            "type": "nominal",
                            "scale": {
                                "range": ['#f77189', '#bb9832', '#50b131', '#36ada4', '#3ba3ec', '#e866f4']
                            }
                        },
                        "opacity": {
                            "condition": {"param": "toggle-curves", "value": 1},
                            "value": 0.2
                        }
                    },
                },
            ];
        if ('values' in m.heaps_curve) {
            layer_values.push({
                "data": m.heaps_curve,
                "mark": {
                    "type": "line",
                    "color": "red"
                },
                "encoding": {
                    "x": x_encoding,
                    "y": {"field": "value", "scale": { "type": "linear" } }
                }
            });
        }
        let yourVlSpec = {
            $schema: 'https://vega.github.io/schema/vega-lite/v6.json',
            "description": 'MultiBar',
            "width": 1000,
            "autosize": {
                "type": "fit",
                "contains": "padding"
            },
            "height": 350,
            "data": m.data,
            "layer": layer_values,
            "resolve": {
                "x": "shared",
                "y": "shared"
            }
        };
        if (!isNaN(m.alpha)) {
            yourVlSpec["title"] = "Growth plot with α=" + m.alpha.toFixed(3);
        }

        function render(scaleType, thisId, vlSpec, add_listeners) {
            const copied_spec = JSON.parse(JSON.stringify(vlSpec)); // deep copy
            if (scaleType == "log") {
                copied_spec.layer[0].encoding.y.scale = { type: "log", domainMin: 1 }; // set scale type
                copied_spec.layer[0].encoding.y2 = { datum: 1 }; // set scale type
            } else {
                copied_spec.layer[0].encoding.y.scale = { type: "linear" }; // set scale type
                if ('y2' in copied_spec.layer[0].encoding) {
                    delete copied_spec.layer[0].encoding[y2]; // set scale type
                }
            }
            let opt = {
                "actions": false,
            };
            vegaEmbed(`#${CSS.escape(thisId)}`, copied_spec, opt).then(({ view, spec, vgSpec }) => {
                if (add_listeners) {
                    // Export PNG
                    let png_button = document.getElementById('btn-download-plot-png-' + m.id);
                    png_button.addEventListener('click', () => {
                        view.toImageURL('png').then(url => {
                            const a = document.createElement('a');
                            a.href = url;
                            a.download = 'visualization.png';
                            a.click();
                        });
                    });

                    // Export SVG
                    let svg_button = document.getElementById('btn-download-plot-svg-' + m.id);
                    svg_button.removeEventListener('click', svg_button);
                    svg_button.addEventListener('click', function svg_button() {
                        view.toImageURL('svg').then(url => {
                            const a = document.createElement('a');
                            a.href = url;
                            a.download = 'visualization.svg';
                            a.click();
                        });
                    });

                    // Open in Vega Editor
                    let vega_editor_button = document.getElementById('btn-download-plot-vega-editor-' + m.id);
                    vega_editor_button.addEventListener('click', () => {
                        post_to_vega_editor(window, {
                            mode: 'vega-lite',
                            spec: JSON.stringify(spec, null, 2),
                            renderer: undefined,
                            config: undefined,
                        });
                    });
                }
            });
        }

        render("linear", id, yourVlSpec, true);
    } else if (element instanceof Line) {
        let l = element;
        let thisId = 'chart-line-' + l.id;
        let mySpec = {
            "$schema": "https://vega.github.io/schema/vega-lite/v6.json",
            "description": "Line",
            "data": l.data,
            "width": 1000,
            "height": 400,
            "layer": [
                {
                    "mark": {
                        "type": "line",
                        "point": {
                            "filled": false,
                            "fill": "white"
                        },
                        "tooltip": true,
                    },
                    "encoding": {
                        "x": {"field": "x", "type": "quantitative", "title": l.x_label},
                        "y": {"field": "y", "type": "quantitative", "title": l.y_label},
                    }
                }
            ]
        };

        if (l.log_x) {
            mySpec.layer[0].encoding.x.scale = { type: "log", nice: false }; // set scale type
            // mySpec.layer[1].encoding.x.scale = { type: "log", nice: false }; // set scale type
            if (!("transform" in mySpec)) {
                mySpec.transform = [];
            }
            mySpec.transform.push({"filter": "datum.x > 0"});
        }
        if (l.log_y) {
            mySpec.layer[0].encoding.y.scale = { type: "log", nice: false }; // set scale type
            // mySpec.layer[1].encoding.y.scale = { type: "log", nice: false }; // set scale type
            if (!("transform" in mySpec)) {
                mySpec.transform = [];
            }
            mySpec.transform.push({"filter": "datum.y > 0"});
        }
        let opt = {
            "actions": false,
        };
        vegaEmbed(`#${CSS.escape(thisId)}`, mySpec, opt).then(({ view, spec, vgSpec }) => {
            // Export PNG
            let png_button = document.getElementById('btn-download-plot-png-' + l.id);
            png_button.addEventListener('click', () => {
                view.toImageURL('png').then(url => {
                    const a = document.createElement('a');
                    a.href = url;
                    a.download = 'visualization.png';
                    a.click();
                });
            });

            // Export SVG
            let svg_button = document.getElementById('btn-download-plot-svg-' + l.id);
            svg_button.removeEventListener('click', svg_button);
            svg_button.addEventListener('click', function svg_button() {
                view.toImageURL('svg').then(url => {
                    const a = document.createElement('a');
                    a.href = url;
                    a.download = 'visualization.svg';
                    a.click();
                });
            });

            // Open in Vega Editor
            let vega_editor_button = document.getElementById('btn-download-plot-vega-editor-' + l.id);
            vega_editor_button.addEventListener('click', () => {
                post_to_vega_editor(window, {
                    mode: 'vega-lite',
                    spec: JSON.stringify(spec, null, 2),
                    renderer: undefined,
                    config: undefined,
                });
            });
        });
    } else if (element instanceof Chromosomal) {
        let c = element;
        let thisId = 'chart-chromosomal-' + c.id;
        // buildPlotDownload(myChart, h.id, fname);
        let mySpec = {
            "$schema": "https://vega.github.io/schema/vega-lite/v6.json",
            "description": "Chromosomal plot",
            "data": c.data,
            "width": 1100,
            "params": [
                {
                    "name": "grid",
                    "select": "interval",
                    "bind": "scales"
                },
                {
                  "name": "metric_select",
                  "value": c.labels[0],
                  "bind": {
                    "input": "select",
                    "options": c.labels,
                    "name": "Variable: "
                  }
                },
                {
                  "name": "color_range",
                  "expr": "metric_select == 'Growth' ? 'redyellowblue' : 'yelloworangered'"
                }
            ],
            "transform": [
                {
                  "fold": c.labels,
                  "as": ["metric_type", "metric_value"]
                },
                {
                  "filter": "datum.metric_type == metric_select"
                }
            ],
            "mark": {
                "type": "rect",
                "clip": true,
                "tooltip": true
            },
            "title": {
                "text": c.sequence,
                "anchor": "start"
            },
            "encoding": {
                "x": {
                    "field": "x",
                    "type": "quantitative",
                    "scale": {"nice": false, "zero": false},
                    "axis": {"domain": false, "grid": false, "ticks": true, "labels": true},
                    "title": "Position in bps"
                },
                "x2": {
                    "field": "x2",
                },
                "color": {
                    "field": "metric_value",
                    "type": "quantitative",
                    "scale": {
                      "scheme": {"expr": "color_range"} ,
                      "type": "linear"
                    }
                }
            },
            "config": {
                "axis": {"grid": true, "tickBand": "extent"},
                "view": {"stroke": "black", "strokeWidth": 2, "cornerRadius": 20}
            }
          };

        function render(scaleType, thisId, vlSpec, addListeners) {
            const copied_spec = JSON.parse(JSON.stringify(vlSpec)); // deep copy
            console.log("Setting scale type: " + scaleType);
            if (scaleType == "log") {
                copied_spec.encoding.color.scale.type = "symlog"; // set scale type
            } else {
                copied_spec.encoding.color.scale.type = "linear"; // set scale type
            }
            let opt = {
                "actions": false,
            };
        vegaEmbed(`#${CSS.escape(thisId)}`, copied_spec, opt).then(({ view, spec, vgSpec }) => {

            if (addListeners) {

            // Export PNG
            let png_button = document.getElementById('btn-download-plot-png-' + c.id);
            png_button.addEventListener('click', () => {
                view.toImageURL('png').then(url => {
                    const a = document.createElement('a');
                    a.href = url;
                    a.download = 'visualization.png';
                    a.click();
                });
            });

            // Export SVG
            let svg_button = document.getElementById('btn-download-plot-svg-' + c.id);
            svg_button.removeEventListener('click', svg_button);
            svg_button.addEventListener('click', function svg_button() {
                view.toImageURL('svg').then(url => {
                    const a = document.createElement('a');
                    a.href = url;
                    a.download = 'visualization.svg';
                    a.click();
                });
            });

            // Open in Vega Editor
            let vega_editor_button = document.getElementById('btn-download-plot-vega-editor-' + c.id);
            vega_editor_button.addEventListener('click', () => {
                post_to_vega_editor(window, {
                    mode: 'vega-lite',
                    spec: JSON.stringify(spec, null, 2),
                    renderer: undefined,
                    config: undefined,
                });
            });
            }
        });
        }

        console.log('btn-logscale-plot-chrom-' + c.id);
        document.getElementById('btn-logscale-plot-chrom-' + c.id).addEventListener('change', (event) => {
            if (event.currentTarget.checked) {
                render("log", thisId, mySpec, false);
            } else {
                render("linear", thisId, mySpec, false);
            }
        });

        if (document.getElementById('btn-logscale-plot-chrom-' + c.id).checked) {
            render("log", thisId, mySpec, true);
        } else {
            render("linear", thisId, mySpec, true);
        }
    } else if (element instanceof Hexbin) {
        let h = element;
        let thisId = 'chart-hexbin-' + h.id;
        // buildPlotDownload(myChart, h.id, fname);
        let mySpec = {
            "$schema": "https://vega.github.io/schema/vega-lite/v6.json",
            "description": "Hexbin",
            "data": h.bins,
            "width": 795,
            "height": 805,
            "params": [{"name": "highlight", "select": "point"}],
            "mark": {
                "type": "text",
                "text": "⬢",
                "size": 81,
                "clip": true,
                "tooltip": true,
            },
            "encoding": {
                "y": {
                    "field": "length",
                    "title": "log10 length in bp",
                    "type": "quantitative",
                    "scale": {"nice": false, "zero": false },
                },
                "x": {
                    "field": "coverage",
                    "type": "quantitative",
                    "scale": {"nice": false, "zero": false },
                },
                "color": {
                    "field": "size",
                    "type": "quantitative",
                    "scale": {"type": "log", "scheme": "bluepurple"}
                },
                "stroke": {
                    "condition": {
                        "param": "highlight",
                        "empty": false,
                        "value": "black"
                    },
                    "value": null
                },
                "opacity": {
                    "condition": {"param": "highlight", "value": 1},
                    "value": 0.5
                },
                "order": {"condition": {"param": "highlight", "value": 1}, "value": 0}
            }
        };

        let opt = {
            "actions": false,
        };
        vegaEmbed(`#${CSS.escape(thisId)}`, mySpec, opt).then(({ view, spec, vgSpec }) => {
            let list_button = document.getElementById('btn-download-node-list-' + h.id);
            list_button.addEventListener('click', () => {
                let ids = new Array();
                if ("vlPoint" in view.signal('highlight')) {
                    ids = view.signal('highlight').vlPoint.or.map((x) => x._vgsid_);
                }
                let table = "";
                ids.forEach((id) => {
                    h.bin_content[id - 1].forEach((dataPoint) => {
                        if (dataPoint != -1) {
                            table += dataPoint + "\t" + id + "\n";
                        } else {
                            table += "ThisBinContainsMoreThan1000Entries\t" + id + "\n";
                        }
                    });
                });
                let blob = new Blob([table], {type: 'text/plain'});
                let a = document.createElement('a');
                a.href = URL.createObjectURL(blob);
                a.download = 'hexbin_nodes_table.tsv';
                a.click();
            });

            view.addSignalListener('highlight', (name, value) => {
                if ("vlPoint" in value) {
                    list_button.removeAttribute('disabled');
                } else {
                    list_button.setAttribute('disabled', '');
                }
            });
            // Export PNG
            let png_button = document.getElementById('btn-download-plot-png-' + h.id);
            png_button.addEventListener('click', () => {
                view.toImageURL('png').then(url => {
                    const a = document.createElement('a');
                    a.href = url;
                    a.download = 'visualization.png';
                    a.click();
                });
            });

            // Export SVG
            let svg_button = document.getElementById('btn-download-plot-svg-' + h.id);
            svg_button.removeEventListener('click', svg_button);
            svg_button.addEventListener('click', function svg_button() {
                view.toImageURL('svg').then(url => {
                    const a = document.createElement('a');
                    a.href = url;
                    a.download = 'visualization.svg';
                    a.click();
                });
            });

            // Open in Vega Editor
            let vega_editor_button = document.getElementById('btn-download-plot-vega-editor-' + h.id);
            vega_editor_button.addEventListener('click', () => {
                post_to_vega_editor(window, {
                    mode: 'vega-lite',
                    spec: JSON.stringify(spec, null, 2),
                    renderer: undefined,
                    config: undefined,
                });
            });
        });

    } else if (element instanceof Heatmap) {
        let h = element;
        let thisId = 'chart-heatmap-' + h.id;
        // buildPlotDownload(myChart, h.id, fname);
        let mySpec = {
            "$schema": "https://vega.github.io/schema/vega-lite/v6.json",
            "description": "Heatmap",
            "data": h.data_set,
            "width": 800,
            "height": 800,
            "mark": {
                "type": "rect",
                "tooltip": true,
            },
            "encoding": {
                "y": {
                    "field": "y",
                    "type": "ordinal",
                    "sort": null,
                    "axis": {
                        "labelOverlap": "greedy"
                    }
                },
                "x": {
                    "field": "x",
                    "type": "ordinal",
                    "sort": null,
                    "axis": {
                        "labelOverlap": "greedy"
                    }
                },
                "color": {
                    "field": "value",
                    "type": "quantitative",
                    "scale": {"range": ["darkred", "white"], "interpolate": "cubehelix", "domainMax": 1.0}
                },
            }
        };

        let opt = {
            "actions": false,
        };
        vegaEmbed(`#${CSS.escape(thisId)}`, mySpec, opt).then(({ view, spec, vgSpec }) => {
            // Export PNG
            let png_button = document.getElementById('btn-download-plot-png-' + h.id);
            png_button.addEventListener('click', () => {
                view.toImageURL('png').then(url => {
                    const a = document.createElement('a');
                    a.href = url;
                    a.download = 'visualization.png';
                    a.click();
                });
            });

            // Export SVG
            let svg_button = document.getElementById('btn-download-plot-svg-' + h.id);
            svg_button.removeEventListener('click', svg_button);
            svg_button.addEventListener('click', function svg_button() {
                view.toImageURL('svg').then(url => {
                    const a = document.createElement('a');
                    a.href = url;
                    a.download = 'visualization.svg';
                    a.click();
                });
            });

            // Open in Vega Editor
            let vega_editor_button = document.getElementById('btn-download-plot-vega-editor-' + h.id);
            vega_editor_button.addEventListener('click', () => {
                post_to_vega_editor(window, {
                    mode: 'vega-lite',
                    spec: JSON.stringify(spec, null, 2),
                    renderer: undefined,
                    config: undefined,
                });
            });
        });
    } else if (element instanceof VegaPlot) {
        let v = element;
        let thisId = 'chart-line-' + v.id;
        let opt = {
            "actions": false,
        };
        vegaEmbed(`#${CSS.escape(thisId)}`, v.jsonContent, opt).then(({ view, spec, vgSpec }) => {
            // Export PNG
            let png_button = document.getElementById('btn-download-plot-png-' + v.id);
            png_button.addEventListener('click', () => {
                view.toImageURL('png').then(url => {
                    const a = document.createElement('a');
                    a.href = url;
                    a.download = 'visualization.png';
                    a.click();
                });
            });

            // Export SVG
            let svg_button = document.getElementById('btn-download-plot-svg-' + v.id);
            svg_button.removeEventListener('click', svg_button);
            svg_button.addEventListener('click', function svg_button() {
                view.toImageURL('svg').then(url => {
                    const a = document.createElement('a');
                    a.href = url;
                    a.download = 'visualization.svg';
                    a.click();
                });
            });

            // Open in Vega Editor
            let vega_editor_button = document.getElementById('btn-download-plot-vega-editor-' + v.id);
            vega_editor_button.addEventListener('click', () => {
                post_to_vega_editor(window, {
                    mode: 'vega-lite',
                    spec: JSON.stringify(spec, null, 2),
                    renderer: undefined,
                    config: undefined,
                });
            });
        });
    } else if (element instanceof DownloadHelper) {
        let d = element;
        document.getElementById('btn-download-plot-' + d.id).addEventListener('click', () => {
            if (d.type == "png") {
                const a = document.createElement('a');
                let png_img = document.getElementById(d.id).getAttribute('src');
                a.href = png_img;
                a.download = 'visualization.png';
                a.click();
            } else if (d.type == "svg") {
                let svgData = document.getElementById(d.id).innerHTML;
                let svgBlob = new Blob([svgData], {type:"image/svg+xml;charset=utf-8"});
                let svgUrl = URL.createObjectURL(svgBlob);
                const a = document.createElement("a");
                a.href = svgUrl;
                a.download = 'visualization.svg';
                a.click();
            }
        });
    }
}

for (let key in objects.tables) {
    let table = objects.tables[key];
    buildTableDownload(table, key, key + '_' + fname);
}
