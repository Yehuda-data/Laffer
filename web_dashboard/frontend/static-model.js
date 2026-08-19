"use strict";

/* Browser implementation of the validated steady-state model for GitHub Pages. */
(function (root) {
  const TOL = 1e-8;
  const BRANCH_JUMP_THRESHOLD = 0.05;
  const KEY_RATES = [0.80, 0.90, 0.95, 0.99];
  const STATUS_INPUT = "INPUT";
  const STATUS_CALIBRATED = "CALIBRATED";
  const STATUS_IMPLIED = "MODEL-IMPLIED";
  const STATUS_OUTPUT = "EQUILIBRIUM OUTPUT";
  const GUESSES = [0.25, 0.05, 0.10, 0.15, 0.20, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65];

  const PRESETS = {presets: [
    {id: "preset_1", label: "Israel — Input assumptions", specification: {name: "Israel — Input assumptions", closure: "s_laffer", calibration: "external", external_balance_convention: "net_imports", kappa_mode: "labor_target", tau_c: 0.18, tau_n: 0.28, tau_k: 0.3, eta: 2, phi: 1, theta: 0.33, delta: 0.02, kappa: 7.301326188544067, n_target: 0.25, k_y: 7.232075502499556, x_y: 0.1637838055265865, R: 1.017941038215124, psi: 1.035, gamma: 1, debt_y: 0.6762598920008294, g_y: 0.2825851135245842, s_y: 0.03086934880871195, external_balance_y: 0.006519397795431243, other_waste_y: 0, grid_min: 0, grid_max: 0.99, grid_step: 0.001}},
    {id: "preset_2", label: "Israel — Model-implied", specification: {name: "Israel — Model-implied", closure: "s_laffer", calibration: "model_implied", external_balance_convention: "net_imports", kappa_mode: "labor_target", tau_c: 0.18, tau_n: 0.28, tau_k: 0.3, eta: 2, phi: 1, theta: 0.14879189287544134, delta: 0.06736487845411664, kappa: 6.730171670285408, n_target: 0.21174431999999999, k_y: 1.6, x_y: 0.1637838055265865, R: 1.017941038215124, psi: 1.035, gamma: 1, debt_y: 0.6762598920008294, g_y: 0.2825851135245842, s_y: 0.08041896050320213, external_balance_y: 0.006519397795431243, other_waste_y: 0, grid_min: 0, grid_max: 0.99, grid_step: 0.001}},
    {id: "preset_3", label: "Paper / US benchmark", specification: {name: "Paper / US benchmark", closure: "s_laffer", calibration: "external", external_balance_convention: "trade_balance", kappa_mode: "kappa", tau_c: 0.04682307692307692, tau_n: 0.2806615384615384, tau_k: 0.3637615384615384, eta: 2, phi: 1, theta: 0.38, delta: 0.07, kappa: 3.4560634481024484, n_target: 0.25000000000000006, k_y: 2.8599488433333895, x_y: 0.2573953959000051, R: 1.04, psi: 1.02, gamma: 1, debt_y: 0.6279230769230769, g_y: 0.17876923076923076, s_y: 0.07616342127299366, external_balance_y: -0.03576923076923077, other_waste_y: 0, grid_min: 0, grid_max: 0.99, grid_step: 0.001}}
  ]};

  const finite = value => Number.isFinite(Number(value));
  const safeRatio = (a, b) => finite(b) && Number(b) !== 0 ? Number(a) / Number(b) : NaN;
  const clean = value => {
    if (Array.isArray(value)) return value.map(clean);
    if (value && typeof value === "object") return Object.fromEntries(Object.entries(value).map(([key, item]) => [key, clean(item)]));
    return typeof value === "number" && !Number.isFinite(value) ? null : value;
  };
  const clone = value => JSON.parse(JSON.stringify(value));

  function validateSpecification(input) {
    const spec = {...input};
    if (spec.calibration === "external" && (spec.theta == null || spec.delta == null)) throw new Error("external calibration requires theta and delta");
    if (spec.calibration === "model_implied" && (spec.k_y == null || spec.x_y == null)) throw new Error("model-implied calibration requires k_y and x_y");
    if (spec.kappa_mode === "kappa" && spec.kappa == null) throw new Error("kappa mode requires kappa");
    if (spec.kappa_mode === "labor_target" && spec.n_target == null) throw new Error("labor-target mode requires n_target");
    if (spec.closure === "s_laffer" && spec.g_y == null) throw new Error("s-Laffer requires baseline g_y");
    if (spec.closure === "g_laffer" && spec.s_y == null) throw new Error("g-Laffer requires baseline s_y");
    if (!(spec.grid_step > 0) || spec.grid_max < spec.grid_min) throw new Error("invalid tax grid");
    return spec;
  }

  function solveRoot(residual) {
    for (const guess of GUESSES) {
      let x = guess;
      let previous = NaN;
      for (let iteration = 0; iteration < 120; iteration += 1) {
        const fx = residual(x);
        if (!finite(fx)) break;
        if (Math.abs(fx) < 1e-10 && x > 0) return x;
        const h = Math.max(1e-7, Math.abs(x) * 1e-6);
        const derivative = (residual(x + h) - residual(Math.max(1e-12, x - h))) / (x + h - Math.max(1e-12, x - h));
        if (!finite(derivative) || Math.abs(derivative) < 1e-14) break;
        let next = x - fx / derivative;
        if (!finite(next) || next <= 0) next = x / 2;
        if (finite(previous) && Math.abs(next - x) > 2 * Math.abs(x - previous)) next = (next + x) / 2;
        previous = x;
        x = next;
      }
      if (x > 0 && finite(residual(x)) && Math.abs(residual(x)) < 1e-9) return x;
    }
    return NaN;
  }

  function steadyS(args) {
    const ky = 1 / ((args.R - 1) / ((1 - args.tau_k) * args.theta) + args.delta / args.theta);
    const yn = Math.pow(args.gamma * Math.pow(ky, args.theta), 1 / (1 - args.theta));
    const xy = (args.psi - 1 + args.delta) * ky;
    const fixedUses = args.g + args.tb + args.waste;
    const alpha = (1 - args.theta) * ((1 - args.tau_n) / (1 + args.tau_c)) * (args.phi / (1 + args.phi));
    const exponent = -(1 + 1 / args.phi);
    const residual = n => {
      const positiveN = n <= 0 ? 1e-12 : n;
      const cy = 1 - xy - fixedUses / yn / positiveN;
      return Math.pow(positiveN, exponent) - args.kappa * (args.eta * cy / alpha + 1 - args.eta);
    };
    const n = solveRoot(residual);
    const y = yn * n;
    const cy = 1 - xy - fixedUses / yn / n;
    const taxY = args.tau_c * cy + args.tau_n * (1 - args.theta) + args.tau_k * (args.theta - args.delta * ky);
    const transfersY = taxY - args.b / y * (args.R - args.psi) - args.g / y - args.waste / y;
    return {y_bar: y, sy_bar: transfersY, ky_bar: ky, xy_bar: xy, cy_bar: cy, n_bar: n, govcons_bar: args.g,
      taxrev_bar: taxY * y, constaxrev_bar: args.tau_c * cy * y, labtaxrev_bar: args.tau_n * (1 - args.theta) * y,
      captaxrev_bar: args.tau_k * (args.theta - args.delta * ky) * y};
  }

  function steadyG(args) {
    const ky = 1 / ((args.R - 1) / ((1 - args.tau_k) * args.theta) + args.delta / args.theta);
    const yn = Math.pow(args.gamma * Math.pow(ky, args.theta), 1 / (1 - args.theta));
    const xy = (args.psi - 1 + args.delta) * ky;
    const incomeTaxY = args.tau_n * (1 - args.theta) + args.tau_k * (args.theta - args.delta * ky);
    const unaffected = args.b * (args.R - args.psi) + args.s - args.tb;
    const consumptionY = n => (1 - xy - incomeTaxY + unaffected / (yn * n)) / (1 + args.tau_c);
    const alpha = (1 - args.theta) * ((1 - args.tau_n) / (1 + args.tau_c)) * (args.phi / (1 + args.phi));
    const exponent = -(1 + 1 / args.phi);
    const residual = n => {
      const positiveN = n <= 0 ? 1e-12 : n;
      return Math.pow(positiveN, exponent) - args.kappa * (args.eta * consumptionY(positiveN) / alpha + 1 - args.eta);
    };
    const n = solveRoot(residual);
    const y = yn * n;
    const cy = consumptionY(n);
    const taxY = args.tau_c * cy + incomeTaxY;
    const tax = taxY * y;
    const g = tax - args.b * (args.R - args.psi) - args.s - args.waste;
    return {y_bar: y, sy_bar: args.s / y, ky_bar: ky, xy_bar: xy, cy_bar: cy, n_bar: n, govcons_bar: g,
      taxrev_bar: tax, constaxrev_bar: args.tau_c * cy * y, labtaxrev_bar: args.tau_n * (1 - args.theta) * y,
      captaxrev_bar: args.tau_k * (args.theta - args.delta * ky) * y};
  }

  function tradeBalanceY(spec) {
    return spec.external_balance_convention === "net_imports" ? -spec.external_balance_y : spec.external_balance_y;
  }

  function structuralValues(spec) {
    if (spec.calibration === "model_implied") {
      const delta = spec.x_y / spec.k_y - (spec.psi - 1);
      const theta = spec.k_y * ((spec.R - 1) / (1 - spec.tau_k) + delta);
      return [theta, delta];
    }
    return [spec.theta, spec.delta];
  }

  function cfeKappa(theta, tauN, tauC, phi, eta, cy, n) {
    const alpha = (1 - theta) * ((1 - tauN) / (1 + tauC)) * (phi / (1 + phi));
    const denominator = eta * cy / alpha + 1 - eta;
    if (!(denominator > 0)) throw new Error(`Invalid CFE kappa denominator: ${denominator}`);
    return Math.pow(n, -(1 + 1 / phi)) / denominator;
  }

  function calibrate(specification) {
    const spec = validateSpecification(specification);
    const [theta, delta] = structuralValues(spec);
    const ky = 1 / ((spec.R - 1) / ((1 - spec.tau_k) * theta) + delta / theta);
    const xy = (spec.psi - 1 + delta) * ky;
    const tbY = tradeBalanceY(spec);
    const incomeTaxY = spec.tau_n * (1 - theta) + spec.tau_k * (theta - delta * ky);
    const cy = spec.closure === "s_laffer"
      ? 1 - xy - spec.g_y - tbY - spec.other_waste_y
      : (1 - xy - incomeTaxY + spec.debt_y * (spec.R - spec.psi) + spec.s_y - tbY) / (1 + spec.tau_c);
    let kappa;
    let n;
    let kappaStatus;
    if (spec.kappa_mode === "labor_target") {
      n = spec.n_target;
      kappa = cfeKappa(theta, spec.tau_n, spec.tau_c, spec.phi, spec.eta, cy, n);
      kappaStatus = STATUS_CALIBRATED;
    } else {
      kappa = spec.kappa;
      const alpha = (1 - theta) * ((1 - spec.tau_n) / (1 + spec.tau_c)) * (spec.phi / (1 + spec.phi));
      const base = kappa * (spec.eta * cy / alpha + 1 - spec.eta);
      n = base > 0 ? Math.pow(base, -spec.phi / (spec.phi + 1)) : NaN;
      kappaStatus = STATUS_INPUT;
    }
    const yn = Math.pow(spec.gamma * Math.pow(ky, theta), 1 / (1 - theta));
    const y = yn * n;
    return {spec, theta, delta, kappa, kappaStatus, n_anchor: n, k_y_anchor: ky, x_y_anchor: xy, c_y_anchor: cy,
      y_anchor: y, b: spec.debt_y * y, g: spec.closure === "s_laffer" ? spec.g_y * y : null,
      s: spec.closure === "g_laffer" ? spec.s_y * y : null, tb: tbY * y, waste: spec.other_waste_y * y};
  }

  function engineSolution(context, tauN = null, tauK = null) {
    const spec = context.spec;
    const args = {b: context.b, tb: context.tb, theta: context.theta, delta: context.delta, kappa: context.kappa,
      waste: context.waste, tau_n: tauN == null ? spec.tau_n : tauN, tau_c: spec.tau_c,
      tau_k: tauK == null ? spec.tau_k : tauK, psi: spec.psi, phi: spec.phi, eta: spec.eta, R: spec.R, gamma: spec.gamma};
    if (spec.closure === "s_laffer") {
      const solution = steadyS({...args, g: context.g});
      return [solution, context.g, solution.sy_bar * solution.y_bar];
    }
    const solution = steadyG({...args, s: context.s});
    return [solution, solution.govcons_bar, context.s];
  }

  function point(context, tauN = null, tauK = null) {
    const spec = context.spec;
    const actualTauN = tauN == null ? spec.tau_n : tauN;
    const actualTauK = tauK == null ? spec.tau_k : tauK;
    try {
      const [solution, g, s] = engineSolution(context, actualTauN, actualTauK);
      const y = solution.y_bar;
      const n = solution.n_bar;
      const ky = solution.ky_bar;
      const xy = solution.xy_bar;
      const cy = solution.cy_bar;
      const k = ky * y;
      const x = xy * y;
      const c = cy * y;
      const w = (1 - context.theta) * y / n;
      const d = context.theta / ky;
      const laborBase = (1 - context.theta) * y;
      const capitalBase = (context.theta - context.delta * ky) * y;
      const total = solution.taxrev_bar;
      const resourceResidual = y - c - x - g - context.tb - context.waste;
      const budgetResidual = total - context.b * (spec.R - spec.psi) - g - s - context.waste;
      const values = [y, n, ky, xy, cy, k, x, c, g, s, total];
      const reasons = [];
      if (!values.every(finite)) reasons.push("solver failure or non-finite equilibrium value");
      if (values.every(finite) && !(n > 0 && n <= 1)) reasons.push("labor is outside 0 < n <= 1");
      if (values.every(finite) && c < 0) reasons.push("consumption is negative");
      if (values.every(finite) && g < 0) reasons.push("government spending is negative");
      if (values.every(finite) && ky < 0) reasons.push("capital/output is negative");
      if (finite(resourceResidual) && Math.abs(resourceResidual) > TOL) reasons.push("resource residual exceeds tolerance");
      if (finite(budgetResidual) && Math.abs(budgetResidual) > TOL) reasons.push("government-budget residual exceeds tolerance");
      const valid = reasons.length === 0;
      return {tau_n: actualTauN, tau_k: actualTauK, point_valid: valid, valid, invalid_reasons: reasons,
        n, y, k, k_y: ky, x, x_y: xy, c, c_y: cy, g, g_y: safeRatio(g, y), s, s_y: safeRatio(s, y),
        m_y: -safeRatio(context.tb, y), tb_y: safeRatio(context.tb, y), w, d,
        labor_tax_base: laborBase, capital_tax_base: capitalBase, consumption_tax_base: c,
        labor_tax_base_y: safeRatio(laborBase, y), capital_tax_base_y: safeRatio(capitalBase, y), consumption_tax_base_y: safeRatio(c, y),
        T_n: solution.labtaxrev_bar, T_k: solution.captaxrev_bar, T_c: solution.constaxrev_bar, T_total: total,
        T_n_y: safeRatio(solution.labtaxrev_bar, y), T_k_y: safeRatio(solution.captaxrev_bar, y),
        T_c_y: safeRatio(solution.constaxrev_bar, y), T_total_y: safeRatio(total, y), resource_residual: resourceResidual,
        government_budget_residual: budgetResidual};
    } catch (error) {
      return {tau_n: actualTauN, tau_k: actualTauK, point_valid: false, valid: false,
        invalid_reasons: [`solver failure: ${error.name}: ${error.message}`]};
    }
  }

  const diagnostic = (level, code, message) => ({level, code, message});
  function baselineDiagnostics(context, item) {
    const spec = context.spec;
    const diagnostics = [];
    const ky = item.k_y ?? context.k_y_anchor;
    const xy = item.x_y ?? context.x_y_anchor;
    if (finite(ky) && ky > 4) diagnostics.push(diagnostic("WARNING", "HIGH_KY", `k/y = ${Number(ky).toPrecision(4)} is unusually high. It follows from the supplied theta, delta, R, and tau_k through equation (15).`));
    if (finite(xy) && xy > 0.30) diagnostics.push(diagnostic("WARNING", "HIGH_XY", `x/y = ${Number(xy).toPrecision(4)} exceeds 0.30.`));
    if (context.delta <= 0) diagnostics.push(diagnostic("WARNING", "NONPOSITIVE_DELTA", `delta = ${context.delta.toPrecision(4)} is non-positive.`));
    else if (context.delta > 0.20) diagnostics.push(diagnostic("WARNING", "HIGH_DELTA", `delta = ${context.delta.toPrecision(4)} exceeds 0.20.`));
    if (!(context.theta > 0 && context.theta < 1)) diagnostics.push(diagnostic("INVALID EQUILIBRIUM", "INVALID_THETA", `theta = ${context.theta.toPrecision(4)} is outside (0, 1).`));
    if (finite(item.n) && !(item.n > 0 && item.n <= 1)) diagnostics.push(diagnostic("INVALID EQUILIBRIUM", "INVALID_LABOR", `n = ${item.n.toPrecision(4)} is outside (0, 1].`));
    if (finite(item.c) && item.c < 0) diagnostics.push(diagnostic("INVALID EQUILIBRIUM", "NEGATIVE_CONSUMPTION", `c = ${item.c.toPrecision(4)} is negative.`));
    if (finite(item.g) && item.g < 0) diagnostics.push(diagnostic("INVALID EQUILIBRIUM", "NEGATIVE_GOVERNMENT", `g = ${item.g.toPrecision(4)} is negative.`));
    if (finite(item.c_y) && Math.abs(item.c_y) < 0.01) diagnostics.push(diagnostic("WARNING", "NEAR_ZERO_CY", `c/y = ${item.c_y.toPrecision(4)} is close to zero.`));
    if (spec.R < spec.psi) diagnostics.push(diagnostic("WARNING", "R_BELOW_PSI", `R = ${spec.R.toPrecision(4)} is below psi = ${spec.psi.toPrecision(4)}.`));
    const beta = Math.pow(spec.psi, spec.eta) / spec.R;
    if (finite(beta) && beta > 1) diagnostics.push(diagnostic("WARNING", "BETA_ABOVE_ONE", `Implied beta = psi^eta/R = ${beta.toPrecision(4)} exceeds 1.`));
    if (!item.valid) diagnostics.push(diagnostic("INVALID EQUILIBRIUM", "SOLVER_OR_ADMISSIBILITY", (item.invalid_reasons || ["Equilibrium is invalid."]).join("; ")));
    if (!diagnostics.length) diagnostics.push(diagnostic("INFO", "ADMISSIBLE", "No configured diagnostic threshold was triggered."));
    return diagnostics;
  }

  function parameterPayload(context) {
    const spec = context.spec;
    const structuralStatus = spec.calibration === "model_implied" ? STATUS_IMPLIED : STATUS_INPUT;
    return {tau_c: {value: spec.tau_c, status: STATUS_INPUT}, tau_n: {value: spec.tau_n, status: STATUS_INPUT},
      tau_k: {value: spec.tau_k, status: STATUS_INPUT}, theta: {value: context.theta, status: structuralStatus},
      delta: {value: context.delta, status: structuralStatus}, eta: {value: spec.eta, status: STATUS_INPUT},
      phi: {value: spec.phi, status: STATUS_INPUT}, kappa: {value: context.kappa, status: context.kappaStatus},
      R: {value: spec.R, status: STATUS_INPUT}, psi: {value: spec.psi, status: STATUS_INPUT}, gamma: {value: spec.gamma, status: STATUS_INPUT}};
  }

  function baseline(specification) {
    const spec = validateSpecification(specification);
    try {
      const context = calibrate(spec);
      const item = point(context);
      const statuses = Object.fromEntries(Object.keys(item).filter(key => !["invalid_reasons", "valid", "point_valid"].includes(key)).map(key => [key, STATUS_OUTPUT]));
      return clean({inputs: spec, parameters: parameterPayload(context), equilibrium: item, statuses,
        diagnostics: baselineDiagnostics(context, item), validity: {valid: Boolean(item.valid), reasons: item.invalid_reasons || []}});
    } catch (error) {
      return clean({inputs: spec, parameters: {}, equilibrium: {}, statuses: {}, diagnostics: [diagnostic("INVALID EQUILIBRIUM", "CALIBRATION_FAILURE", `${error.name}: ${error.message}`)], validity: {valid: false, reasons: [error.message]}});
    }
  }

  function grid(spec, baselineRate) {
    const count = Math.round((spec.grid_max - spec.grid_min) / spec.grid_step);
    const values = [];
    for (let index = 0; index <= count; index += 1) {
      const value = spec.grid_min + index * spec.grid_step;
      if (value >= spec.grid_min - 1e-12 && value <= spec.grid_max + 1e-12) values.push(Number(value.toFixed(10)));
    }
    if (!values.some(value => Math.abs(value - baselineRate) <= 1e-12)) values.push(baselineRate);
    return [...new Set(values)].sort((a, b) => a - b);
  }

  function curve(specification, kind) {
    const spec = validateSpecification(specification);
    if (kind === "capital" && spec.closure === "g_laffer") throw new Error("Capital-tax g-Laffer not yet implemented.");
    const context = calibrate(spec);
    const base = point(context);
    const baselineRate = kind === "labor" ? spec.tau_n : spec.tau_k;
    const points = grid(spec, baselineRate).map(rate => kind === "labor" ? point(context, rate, null) : point(context, null, rate));
    const baselineIndex = points.findIndex(item => Math.abs((kind === "labor" ? item.tau_n : item.tau_k) - baselineRate) <= 1e-10);
    if (kind === "labor" && spec.closure === "g_laffer") {
      const connected = Array(points.length).fill(false);
      if (baselineIndex >= 0 && points[baselineIndex].point_valid) {
        connected[baselineIndex] = true;
        for (let index = baselineIndex + 1; index < points.length; index += 1) connected[index] = connected[index - 1] && Boolean(points[index].point_valid);
        for (let index = baselineIndex - 1; index >= 0; index -= 1) connected[index] = connected[index + 1] && Boolean(points[index].point_valid);
      }
      points.forEach((item, index) => { item.valid = connected[index]; if (item.point_valid && !item.valid) item.invalid_reasons.push("valid root is not connected to the baseline branch"); });
    }
    for (const item of points) {
      for (const key of ["n", "y", "k", "c"]) item[`${key}_index`] = safeRatio(item[key], base[key]) * 100;
      for (const key of ["T_n", "T_k", "T_c", "T_total"]) item[`${key}_index`] = safeRatio(item[key], base.T_total) * 100;
    }
    const validPoints = points.filter(item => item.valid && finite(item.T_total_index));
    const peak = validPoints.reduce((best, item) => !best || item.T_total_index > best.T_total_index ? item : best, null);
    const diagnostics = baselineDiagnostics(context, base);
    const rates = validPoints.map(item => kind === "labor" ? item.tau_n : item.tau_k);
    const summary = {baseline_tax: baselineRate, peak_tax: peak ? (kind === "labor" ? peak.tau_n : peak.tau_k) : null,
      peak_revenue: peak?.T_total_index ?? null};
    if (kind === "labor") {
      const jumps = [];
      const failures = [];
      for (let index = 1; index < points.length; index += 1) {
        const previous = points[index - 1];
        const current = points[index];
        if (!finite(current.n)) failures.push(current.tau_n);
        if (finite(previous.n) && finite(current.n) && Math.abs(current.n - previous.n) > BRANCH_JUMP_THRESHOLD) jumps.push({tau_n: current.tau_n, delta_n: Math.abs(current.n - previous.n)});
      }
      const firstRate = predicate => points.find(predicate)?.tau_n ?? null;
      Object.assign(summary, {valid_tau_min: rates.length ? Math.min(...rates) : null, valid_tau_max: rates.length ? Math.max(...rates) : null,
        first_g_nonpositive: firstRate(item => finite(item.g) && item.g <= 0), first_n_at_least_one: firstRate(item => finite(item.n) && item.n >= 1),
        first_c_nonpositive: firstRate(item => finite(item.c) && item.c <= 0), solver_failure_rates: failures, possible_branch_switches: jumps});
      if (jumps.length) diagnostics.push(diagnostic("WARNING", "POSSIBLE_BRANCH_SWITCH", "Possible branch/root switch: adjacent labor jump exceeds 0.05."));
      const disconnected = points.filter(item => item.point_valid && !item.valid).length;
      if (disconnected) diagnostics.push(diagnostic("WARNING", "DISCONNECTED_ROOTS", `${disconnected} raw valid roots are outside the baseline-connected branch.`));
    }
    return clean({inputs: spec, parameters: parameterPayload(context), baseline: base, curve: points, summary, diagnostics,
      validity: {valid: Boolean(validPoints.length), valid_points: validPoints.length, total_points: points.length}});
  }

  function sensitivity(request) {
    const allowed = new Set(["theta", "delta", "R", "psi", "eta", "phi", "kappa", "k_y", "x_y", "g_y", "s_y", "m_y"]);
    if (!allowed.has(request.parameter)) throw new Error(`unsupported sensitivity parameter: ${request.parameter}`);
    const spec = validateSpecification(request.specification);
    if (spec.calibration === "model_implied" && ["theta", "delta"].includes(request.parameter)) throw new Error("theta and delta are model-implied and cannot be edited directly");
    if (request.parameter === "g_y" && spec.closure !== "s_laffer") throw new Error("g_y sensitivity is available only for s-Laffer");
    if (request.parameter === "s_y" && spec.closure !== "g_laffer") throw new Error("s_y sensitivity is available only for g-Laffer");
    const scenarios = [];
    for (let index = 0; index < request.scenarios; index += 1) {
      const value = request.minimum + (request.maximum - request.minimum) * index / (request.scenarios - 1);
      const updated = {...spec};
      if (request.parameter === "m_y") { updated.external_balance_convention = "net_imports"; updated.external_balance_y = value; }
      else updated[request.parameter] = value;
      if (request.parameter === "kappa") { updated.kappa_mode = "kappa"; updated.kappa = value; }
      scenarios.push({label: `${request.parameter} = ${value.toPrecision(6)}`, value, result: curve(updated, "labor")});
    }
    return clean({parameter: request.parameter, scenarios});
  }

  function compare(request) {
    const scenario = spec => {
      const baselineResult = baseline(spec);
      const curveResult = curve(spec, "labor");
      const keyRates = KEY_RATES.map(rate => {
        const row = curveResult.curve.reduce((best, item) => !best || Math.abs(item.tau_n - rate) < Math.abs(best.tau_n - rate) ? item : best, null) || {};
        return {tau_n: rate, n: row.valid ? row.n : null, T_total_index: row.valid ? row.T_total_index : null, valid: Boolean(row.valid), reasons: row.invalid_reasons || []};
      });
      return {baseline: baselineResult, curve: curveResult, key_rates: keyRates};
    };
    return clean({scenario_a: scenario(request.scenario_a), scenario_b: scenario(request.scenario_b)});
  }

  function equations(closure) {
    const items = [
      {name: "Capital/output", latex: "k/y=[(R-1)/(theta(1-tau_k))+delta/theta]^{-1}", source: "laffer_model.steady / equation (15)"},
      {name: "Capital accumulation", latex: "x/y=(psi-1+delta)k/y", source: "laffer_model.steady"},
      {name: "Production", latex: "y/n=[gamma(k/y)^theta]^{1/(1-theta)}", source: "laffer_model.steady"},
      {name: "Tax revenue", latex: "T/y=tau_c(c/y)+tau_n(1-theta)+tau_k[theta-delta(k/y)]", source: "laffer_model.steady and laffer_model_g.steady_g"},
      {name: "CFE labor", latex: "n^{-(1+1/phi)}=kappa[eta(c/y)/alpha+1-eta]", source: "existing CFE residual solvers"}
    ];
    if (closure === "s_laffer") items.push({name: "s-Laffer closure", latex: "s/y=T/y-[b/y](R-psi)-g/y-q/y; g fixed in levels", source: "laffer_model.steady"});
    else items.push({name: "g-Laffer consumption", latex: "c/y=(1+tau_c)^{-1}{1-x/y-tau_n(1-theta)-tau_k[theta-delta(k/y)]+[b(R-psi)+s-tb]/y}", source: "laffer_model_g._cy_g, equation (19)"}, {name: "g-Laffer closure", latex: "g=T-b(R-psi)-s-q; s fixed in levels", source: "laffer_model_g.steady_g / Proposition 3"});
    return {closure, equations: items};
  }

  async function request(path, options = {}) {
    const body = options.body ? JSON.parse(options.body) : null;
    if (path === "/api/presets") return clone(PRESETS);
    if (path === "/api/baseline") return baseline(body);
    if (path === "/api/laffer/labor") return curve(body, "labor");
    if (path === "/api/laffer/capital") return curve(body, "capital");
    if (path === "/api/sensitivity") return sensitivity(body);
    if (path === "/api/compare") return compare(body);
    if (path.startsWith("/api/equations/")) return equations(path.split("/").pop());
    if (path === "/api/health") return {status: "ok"};
    throw new Error(`Unknown static API path: ${path}`);
  }

  const api = {request, baseline, laborCurve: spec => curve(spec, "labor"), capitalCurve: spec => curve(spec, "capital"), sensitivity, compare, equations, presets: () => clone(PRESETS)};
  root.LafferStaticApi = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window !== "undefined" ? window : globalThis);
