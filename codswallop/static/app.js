/* CODSWALLOP -- shared behaviour present on every page.
 *
 * Just the theme toggle so far. The family page's machinery lives in family.js and map.js,
 * because it needs the family data and this does not.
 */
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
