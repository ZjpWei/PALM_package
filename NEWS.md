# PALM 0.1.1

* The count-matrix argument of `palm()` and `palm.null.model()` is renamed from
  `rel.abd` to `feature.counts`. `rel.abd` is still accepted, with a deprecation
  warning, so existing code keeps running; it gives identical results.
* Includes the stability updates made after the 0.1.0 release, so results
  may differ from 0.1.0.
* Fixed the help page of `palm()`: a duplicated paragraph with an unescaped `%`
  truncated it, dropping the value, examples, see-also and author sections.

# PALM 0.1.0

* Version used in Wei et al., *Fast and reliable association discovery in
  large-scale microbiome studies and meta-analyses using PALM*. Install it with
  `devtools::install_github("ZjpWei/PALM_package@Package_v1")`.
