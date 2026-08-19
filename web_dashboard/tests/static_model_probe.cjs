"use strict";

const fs = require("node:fs");
const api = require("../frontend/static-model.js");
const specifications = JSON.parse(fs.readFileSync(0, "utf8"));
const output = specifications.map(specification => {
  const baseline = api.baseline(specification);
  const labor = api.laborCurve(specification);
  const capital = specification.closure === "s_laffer" ? api.capitalCurve(specification) : null;
  return {
    baseline: {
      valid: baseline.validity.valid,
      n: baseline.equilibrium.n,
      y: baseline.equilibrium.y,
      k_y: baseline.equilibrium.k_y,
      c_y: baseline.equilibrium.c_y,
      T_total_y: baseline.equilibrium.T_total_y,
    },
    labor: {
      peak_tax: labor.summary.peak_tax,
      peak_revenue: labor.summary.peak_revenue,
      valid_points: labor.validity.valid_points,
    },
    capital: capital ? {
      peak_tax: capital.summary.peak_tax,
      peak_revenue: capital.summary.peak_revenue,
      valid_points: capital.validity.valid_points,
    } : null,
  };
});
process.stdout.write(JSON.stringify(output));
