/* CODSWALLOP -- shared behaviour present on every page.
 *
 * The theme toggle, and the footer's "what is already filed" line. The family page's
 * machinery lives in family.js and map.js, because it needs the family data and this does not.
 */

/* What the archive already holds, in the footer of every page.
 *
 * This also IS the per-app hit beacon for the mdeller.com launcher, which counts a page view
 * by watching the access log for a request only a rendering browser makes: a scanner fetches
 * the HTML and stops. The beacon in apps.json is `^/api/stats\b`, so this fetch has to exist
 * on every page and has to run after render, or the app's hit count reads zero forever with
 * nothing to indicate why. It was declared as the beacon before anything called it.
 *
 * A static asset would be the easier beacon and a worse one: /static/ is served
 * `max-age=31536000, immutable`, so a returning reader never re-requests it and the count
 * silently undercounts.
 */
(function () {
  "use strict";

  fetch("/api/stats", { headers: { "Accept": "application/json" } })
    .then(function (r) { return r.ok ? r.json() : null; })
    .then(function (s) {
      if (!s || !s.families) return;
      var box = document.getElementById("filedStat");
      if (!box) return;
      var n = function (v) { return Number(v).toLocaleString("en-GB"); };
      box.textContent = "Filed so far: " + n(s.families) + " famil" +
        (s.families === 1 ? "y" : "ies") + ", " + n(s.entries) + " PDB entries, " +
        n(s.entities) + " polymer entities cached.";
      box.hidden = false;
    })
    .catch(function () { /* the footer line is decoration; never break a page for it */ });
})();

(function () {
  "use strict";

  var KEY = "codswallop-theme";
  var root = document.documentElement;
  var btn = document.getElementById("themeToggle");
  if (!btn) return;

  function current() {
    var set = root.getAttribute("data-theme");
    if (set === "light" || set === "dark") return set;
    // Nothing chosen yet: report what the OS preference is actually rendering, so the
    // first click flips what the reader sees rather than what a default says they see.
    return window.matchMedia("(prefers-color-scheme: light)").matches ? "light" : "dark";
  }

  function label() {
    var next = current() === "dark" ? "light" : "dark";
    btn.setAttribute("title", "Switch to the " + next + " archive");
    btn.setAttribute("aria-label", "Switch to the " + next + " archive");
  }

  btn.addEventListener("click", function () {
    var next = current() === "dark" ? "light" : "dark";
    root.setAttribute("data-theme", next);
    try { localStorage.setItem(KEY, next); } catch (e) { /* private browsing */ }
    label();
  });

  label();
})();
