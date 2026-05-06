// =============================================================================
// turn 145 — Server unify, atlas-side: standalone button + popup
// =============================================================================
// SPEC: docs/SERVER_AUDIT_2026-05-05.md §"Status button (top bar, turn 145)"
//
// Verifies (source-level, no JS sandbox needed for most):
//   - The server status chip has been LIFTED OUT of the data dropdown
//   - A standalone top-bar button (#atlasServerStandaloneBtn) exists with
//     the right anatomy: dot, label, hover title, neutral default style
//   - The popup overlay (#atlasServerPopupOverlay) exists, hidden by default
//   - The atlasServer object captures lastHealthBody on probe
//   - The init function is exposed on window
//
// Plus a small JS-sandbox section for the popup's HTML-rendering logic
// to confirm the two snippet contents are exactly what the user copies.
// =============================================================================

const fs = require('fs');
const path = require('path');
const vm = require('vm');

const ATLAS_PATH = path.resolve(__dirname, '..', 'Inversion_atlas.html');
const html = fs.readFileSync(ATLAS_PATH, 'utf8');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
  if (cond) { pass++; console.log('  PASS ' + name); }
  else { fail++; console.log('  FAIL ' + name + (detail ? ' :: ' + detail : '')); }
}

// ============================================================================
// 1. Standalone button DOM
// ============================================================================
console.log('\n=== 1. Standalone button DOM ===');
ok('#atlasServerStandaloneBtn declared',
   /id="atlasServerStandaloneBtn"/.test(html));
ok('button has class atlas-server-standalone',
   /id="atlasServerStandaloneBtn"[^>]*class="atlas-server-standalone"/.test(html));
ok('contains the colored dot span',
   /id="atlasServerStandaloneDot"/.test(html));
ok('contains the label span "server"',
   /id="atlasServerStandaloneLabel"[^>]*>server<\/span>/.test(html));
ok('button has descriptive title',
   /id="atlasServerStandaloneBtn"[^>]*title="Atlas server status\. Click for the start command and config path\."/.test(html));
ok('lives outside the yellow data dropdown',
   (() => {
     const dataPanelEnd = html.indexOf('id="atlasToolsGroup">');
     const dataPanelClose = html.indexOf('</div>\n  </div>', dataPanelEnd);
     const standalone = html.indexOf('id="atlasServerStandaloneBtn"');
     // standalone button must come AFTER the data panel closes
     return standalone > 0 && dataPanelClose > 0 && standalone > dataPanelClose;
   })());

// ============================================================================
// 2. Data dropdown no longer contains the dynamic badge target comment
// ============================================================================
console.log('\n=== 2. Data dropdown cleaned of legacy badge ===');
ok('legacy "appended dynamically by atlasServer" comment is gone',
   !/Server status badge appended dynamically/.test(html));
ok('legacy #atlasServerBadge id no longer in DOM',
   (html.match(/id="atlasServerBadge"/g) || []).length === 0);
ok('#activeSamplesBadge still in data dropdown (sanity — only server chip moved)',
   /id="activeSamplesBadge"/.test(html));

// ============================================================================
// 3. Popup overlay shell
// ============================================================================
console.log('\n=== 3. Popup overlay shell ===');
ok('#atlasServerPopupOverlay declared',
   /id="atlasServerPopupOverlay"/.test(html));
ok('overlay hidden by default (display: none)',
   /id="atlasServerPopupOverlay"[\s\S]{0,300}?display:\s*none/.test(html));
ok('overlay covers viewport (100vw + 100vh)',
   /id="atlasServerPopupOverlay"[\s\S]{0,400}?100vw[\s\S]{0,200}?100vh/.test(html));
ok('overlay has high z-index (>= 9000)',
   /id="atlasServerPopupOverlay"[\s\S]{0,400}?z-index:\s*9\d{3}/.test(html));
ok('inner #atlasServerPopup container',
   /id="atlasServerPopup"/.test(html));
ok('popup has role=dialog',
   /id="atlasServerPopup"[\s\S]{0,200}?role="dialog"/.test(html));
ok('popup has aria-labelledby',
   /id="atlasServerPopup"[\s\S]{0,200}?aria-labelledby="atlasServerPopupTitle"/.test(html));

// ============================================================================
// 4. Init function + window export
// ============================================================================
console.log('\n=== 4. Init function definition + exports ===');
ok('_atlasServerInitStandaloneBtn defined',
   /function _atlasServerInitStandaloneBtn\(\)/.test(html));
ok('legacy _atlasServerInitBadge kept as no-op alias',
   /function _atlasServerInitBadge\(\)\s*\{\s*return _atlasServerInitStandaloneBtn\(\)/.test(html));
ok('window._atlasServerInitStandaloneBtn export',
   /window\._atlasServerInitStandaloneBtn\s*=\s*_atlasServerInitStandaloneBtn/.test(html));
ok('window._atlasServerInitBadge backwards-compat export',
   /window\._atlasServerInitBadge\s*=\s*_atlasServerInitBadge/.test(html));
ok('Init wired on DOMContentLoaded or setTimeout',
   /DOMContentLoaded[\s\S]{0,80}?_atlasServerInitStandaloneBtn|setTimeout\(_atlasServerInitStandaloneBtn/.test(html));

// ============================================================================
// 5. Popup content — start + config snippets present in source
// ============================================================================
console.log('\n=== 5. Popup snippet contents ===');
ok('start snippet is exactly "./run_atlas.sh"',
   /['"]\.\/run_atlas\.sh['"]/.test(html));
ok('cross-platform alt is "python3 run_atlas.py"',
   /python3 run_atlas\.py/.test(html));
ok('config snippet copies the example yaml',
   /cp server_turn1\/popstats_server\.config\.example\.yaml/.test(html));
ok('config snippet target is the canonical name',
   /server_turn1\/popstats_server\.config\.yaml/.test(html));
ok('start snippet section labelled START',
   /start\s*&mdash;\s*paste in a terminal/.test(html));
ok('config snippet section labelled CONFIG',
   /config\s*&mdash;\s*one-time/.test(html));

// ============================================================================
// 6. Popup interactions wired
// ============================================================================
console.log('\n=== 6. Popup interaction surfaces ===');
ok('close button id present',
   /id="atlasServerPopupClose"/.test(html));
ok('recheck button id present',
   /id="atlasServerPopupRecheck"/.test(html));
ok('two copy-buttons (one per snippet) using data-copy-target attribute',
   (html.match(/data-copy-target="atlasServerStartSnippet"|data-copy-target="atlasServerCfgSnippet"/g) || []).length === 2);
ok('Esc key closes popup',
   /e\.key === 'Escape'[\s\S]{0,120}?_closePopup/.test(html));
ok('click outside (overlay click) closes popup',
   /overlay\.addEventListener\('click'[\s\S]{0,200}?_closePopup/.test(html));
ok('clipboard fallback path uses execCommand',
   /document\.execCommand\('copy'\)/.test(html));
ok('navigator.clipboard.writeText preferred path',
   /navigator\.clipboard\.writeText/.test(html));

// ============================================================================
// 7. atlasServer captures lastHealthBody
// ============================================================================
console.log('\n=== 7. lastHealthBody capture ===');
ok('lastHealthBody initialized to null in atlasServer',
   /atlasServer = \{[\s\S]{0,800}?lastHealthBody:\s*null/.test(html));
ok('isAvailable assigns lastHealthBody from response JSON',
   /this\.lastHealthBody\s*=\s*await resp\.json\(\)/.test(html));
ok('lastHealthBody reset on probe failure',
   /catch \(_\) \{\s*this\.lastHealthBody\s*=\s*null/.test(html));
ok('setUrl resets lastHealthBody',
   /setUrl\([\s\S]{0,200}?this\.lastHealthBody\s*=\s*null/.test(html));

// ============================================================================
// 8. Three-state styling: ready / partial / down
// ============================================================================
console.log('\n=== 8. Three-state styling ===');
ok('ready palette uses #3cc08a (green)',
   /ready:\s*\{\s*dot:\s*'#3cc08a'/.test(html));
ok('partial palette uses #e0b042 (amber)',
   /partial:\s*\{\s*dot:\s*'#e0b042'/.test(html));
ok('down palette uses neutral grey',
   /down:\s*\{\s*dot:\s*'#6b7388'/.test(html));
ok('partial state branch reads subsystems.popstats.ready === false',
   /subsystems\.popstats[\s\S]{0,80}?ready === false/.test(html));
ok('partial requires file subsystem ready === true',
   /subsystems\.file[\s\S]{0,80}?ready === true/.test(html));

// ============================================================================
// 9. Sandbox the popup-html builder & verify exact snippet text
// ============================================================================
console.log('\n=== 9. Popup HTML builder — sandboxed render ===');
{
  // Extract the standalone-init function body so we can call _popupHtml
  // in isolation. The function relies on a few outer references
  // (atlasServer, _atlasServerEsc) which we mock.
  const initRe = /function _atlasServerInitStandaloneBtn\(\)\s*\{[\s\S]*?\n\}/;
  const m = html.match(initRe);
  const escRe = /function _atlasServerEsc\(s\)\s*\{[\s\S]*?\n\}/;
  const escM = html.match(escRe);
  if (!m || !escM) {
    ok('extract popup builder', false, 'regex failed');
  } else {
    // Build a sandbox where atlasServer.url is set to a known string,
    // status='available', lastHealthBody is null (so state→'ready').
    const ctx = vm.createContext({
      atlasServer: { url: 'http://127.0.0.1:8765', status: 'available', lastHealthBody: null,
                     isAvailable: async () => true },
      document: {
        getElementById: () => null,   // bail out before DOM-touching
        readyState: 'complete',
        addEventListener: () => {},
        body: { appendChild: () => {}, removeChild: () => {} },
        createElement: () => ({ style: {}, focus: () => {}, click: () => {} }),
      },
      navigator: {},
      console,
      setTimeout: setTimeout,
    });
    // Run escape + init source. Since both `getElementById` returns null,
    // the early return at the top fires — but we still want _popupHtml
    // accessible. So we extract _popupHtml differently: run the init
    // wrapped in a way that exposes it.
    //
    // Actually simpler: parse out the inner _popupHtml fn directly.
    const popupRe = /function _popupHtml\(\)\s*\{[\s\S]*?return\s*\[[\s\S]*?\]\.join\(''\);\s*\}/;
    const popupM = m[0].match(popupRe);
    ok('_popupHtml() function extractable', !!popupM);
    if (popupM) {
      // Build a tiny harness that defines _atlasServerEsc + _popupHtml
      // and calls it.
      const harness =
        escM[0] + '\n' +
        // Outer scope provides the function 'state' look-up. The inner
        // function uses atlasServer.url etc.
        'function _currentState() { return ' + JSON.stringify('ready') + '; }\n' +
        popupM[0] + '\n' +
        '__OUT = _popupHtml();';
      const sandboxCtx = vm.createContext({
        atlasServer: { url: 'http://127.0.0.1:8765', status: 'available', lastHealthBody: null },
        __OUT: null,
      });
      try {
        vm.runInContext(harness, sandboxCtx);
        const out = sandboxCtx.__OUT;
        ok('popup HTML returns string',                typeof out === 'string' && out.length > 200);
        ok('popup HTML contains ./run_atlas.sh',       out.indexOf('./run_atlas.sh') >= 0);
        ok('popup HTML contains python3 run_atlas.py', out.indexOf('python3 run_atlas.py') >= 0);
        ok('popup HTML contains config copy command',
           out.indexOf('cp server_turn1/popstats_server.config.example.yaml') >= 0);
        ok('popup HTML escapes the URL',
           out.indexOf('http://127.0.0.1:8765') >= 0);
        ok('popup HTML has Recheck button',
           /id="atlasServerPopupRecheck"[\s\S]{0,400}?Recheck/.test(out));
        ok('popup HTML has close X',
           /id="atlasServerPopupClose"/.test(out));
        ok('popup HTML has both copy buttons',
           (out.match(/data-copy-target="atlasServerStartSnippet"|data-copy-target="atlasServerCfgSnippet"/g) || []).length === 2);
        ok('popup HTML has dialog title',
           /id="atlasServerPopupTitle"[\s\S]{0,300}?Atlas server/.test(out));
        ok('popup escapes <script> in URL (XSS guard)',
           (() => {
             // Re-run with malicious URL
             const ctx2 = vm.createContext({
               atlasServer: { url: 'http://evil"<script>alert(1)</script>',
                              status: 'available', lastHealthBody: null },
               __OUT: null,
             });
             vm.runInContext(harness, ctx2);
             const o = ctx2.__OUT;
             return o.indexOf('<script>alert(1)') === -1 &&
                    o.indexOf('&lt;script&gt;') >= 0;
           })());
      } catch (e) {
        ok('popup harness ran',          false, e.message);
      }
    }
  }
}

// ============================================================================
// 10. Server merge — server-side surface checks (source-level)
// ============================================================================
console.log('\n=== 10. Merged server source surface ===');
const POPSTATS_PY = path.resolve(__dirname, '..', 'server_turn1', 'popstats_server.py');
let psSrc = '';
try { psSrc = fs.readFileSync(POPSTATS_PY, 'utf8'); } catch (_) {}
ok('popstats_server.py present',                 psSrc.length > 1000);
ok('top-of-file mentions "merged" / turn 145',   /merged|turn 145/i.test(psSrc));
ok('LANTA reference removed from header',
   !/Sits on the cluster \(LANTA\)/.test(psSrc));
ok('PROJECT_ROOT module-level state added',      /PROJECT_ROOT:\s*Optional\[Path\]\s*=\s*None/.test(psSrc));
ok('_bootstrap_file function defined',           /def _bootstrap_file\(project_root: Path\)/.test(psSrc));
ok('_ensure_project_root function defined',      /def _ensure_project_root\(\)/.test(psSrc));
ok('_safe_project_path function defined',        /def _safe_project_path\(rel: str\)/.test(psSrc));
ok('GET /file/{path:path} route registered',     /@app\.get\("\/file\/\{path:path\}"\)/.test(psSrc));
ok('POST /file/{path:path} route registered',    /@app\.post\("\/file\/\{path:path\}"\)/.test(psSrc));
ok('POST /compute/{name} route registered',      /@app\.post\("\/compute\/\{name\}"\)/.test(psSrc));
ok('COMPUTE_REGISTRY exposed',                   /COMPUTE_REGISTRY:\s*Dict\[str, Any\]\s*=\s*\{/.test(psSrc));
ok('echo + list_files computes registered',
   /['"]echo['"]:\s*_compute_echo/.test(psSrc) && /['"]list_files['"]:\s*_compute_list_files/.test(psSrc));
ok('subsystems.popstats.ready field shape',      /"ready":\s*True/.test(psSrc) && /"reason":/.test(psSrc));
ok('main() --config is now optional',            /["']--config["']\s*,\s*default=None/.test(psSrc));
ok('main() accepts --project-root',              /["']--project-root["']/.test(psSrc));
ok('main() sets ATLAS_PROJECT_ROOT env before uvicorn',
   /os\.environ\["ATLAS_PROJECT_ROOT"\]\s*=\s*str\(PROJECT_ROOT\)/.test(psSrc));

// ============================================================================
// 11. Launchers + audit doc + duplicate dir gone
// ============================================================================
console.log('\n=== 11. Launchers + audit doc + cleanup ===');
const ATLAS_DIR = path.resolve(__dirname, '..');
ok('run_atlas.sh exists',
   fs.existsSync(path.join(ATLAS_DIR, 'run_atlas.sh')));
ok('run_atlas.sh starts with #!/usr/bin/env bash',
   fs.readFileSync(path.join(ATLAS_DIR, 'run_atlas.sh'), 'utf8').startsWith('#!/usr/bin/env bash'));
ok('run_atlas.py exists',
   fs.existsSync(path.join(ATLAS_DIR, 'run_atlas.py')));
ok('run_atlas.py starts with #!/usr/bin/env python3',
   fs.readFileSync(path.join(ATLAS_DIR, 'run_atlas.py'), 'utf8').startsWith('#!/usr/bin/env python3'));
ok('audit doc shipped at docs/SERVER_AUDIT_2026-05-05.md',
   fs.existsSync(path.join(ATLAS_DIR, 'docs', 'SERVER_AUDIT_2026-05-05.md')));
ok('byte-duplicate server_turn11c_ld_fast/ deleted',
   !fs.existsSync(path.join(ATLAS_DIR, 'server_turn11c_ld_fast')));
ok('LD endpoint test relocated to server_turn1/test_ld_endpoint.py',
   fs.existsSync(path.join(ATLAS_DIR, 'server_turn1', 'test_ld_endpoint.py')));
ok('atlas_server.py is now a thin shim (< 200 lines)',
   (() => {
     const s = fs.readFileSync(path.join(ATLAS_DIR, 'atlas_server.py'), 'utf8');
     return s.length > 0 && s.indexOf('compatibility shim') >= 0 &&
            s.split('\n').length < 200;
   })());

// ============================================================================
// Final tally
// ============================================================================
console.log('\n=============================================================');
console.log('  ' + pass + ' / ' + (pass + fail) + ' tests passed');
console.log('=============================================================');
process.exit(fail === 0 ? 0 : 1);
