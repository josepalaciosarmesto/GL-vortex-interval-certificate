# Manuscript notes for GL-arxiv-final.tex

This note reflects the revised certificate target after the first CAPD run. The old target 0.44536 should not be reported unless a sharper transcript actually proves it.

## Certified target to report

Use the following target in Section 6:

~~~tex
\operatorname{Tr}(\mathcal H_1-2)_-^2 < 0.4456.
~~~

Consequently,

~~~tex
\lambda_0 > 2-\sqrt{0.4456} > 1.3324.
~~~

Do not report the old bound 0.44536 or the old rounded lower bound 1.3326 from the current certificate.

## Concrete TeX edits

In the macro block around the Lieb-Thirring subsection, change:

~~~tex
\def\bottom{1.3326}
\def\turning{0.614489}
~~~

to:

~~~tex
\def\bottom{1.3324}
\def\turning{0.6140}
~~~

In the compact integral display, replace the old heuristic line

~~~tex
\left. \frac{1}{2(\gamma+1)} A_\gamma \right|_{\gamma = 2} \approx 0.44515
~~~

by an interval statement such as

~~~tex
\left. \frac{1}{2(\gamma+1)} A_\gamma \right|_{\gamma = 2} < 0.4454.
~~~

Because the final estimate is checked by the verifier as a combined interval, avoid deriving the final target only by adding rounded display bounds. Instead replace

~~~tex
\mathrm{Tr}(\mathcal H_1-2)_-^2 \le 0.44536
~~~

by

~~~tex
\mathrm{Tr}(\mathcal H_1-2)_-^2 < 0.4456.
~~~

Finally replace

~~~tex
2-\sqrt{0.44536}\approx 1.3326
~~~

by

~~~tex
2-\sqrt{0.4456}>1.3324.
~~~

## Suggested replacement paragraph

~~~tex
The finite numerical assertions in this subsection were checked with outward
interval arithmetic. Starting from the Frobenius expansion at \(r=0.1\), the
CAPD verifier propagates an interval box for the exact solution of
\[
  U''+\frac1rU'-\frac1{r^2}U+(1-U^2)U=0
\]
for all \(U'(0)\) in the certified shooting interval
\[
  U'(0)\in[0.583189495860,\,0.583189495862].
\]
It verifies the endpoint shooting alternatives, positivity of
\(r^{-2}-3(1-U^2)\) up to \(r_0=0.6140\), and the boundary inequality
\(U(7)\ge 1-0.55\,7^{-2}\). The compact quadrature is checked by interval
arithmetic, and after adding the explicit tail allowance the verifier gives
\[
  \operatorname{Tr}(\mathcal H_1-2)_-^2 < 0.4456.
\]
Therefore
\[
  \lambda_0>2-\sqrt{0.4456}>1.3324.
\]
~~~

## Why no eigenvalue computation is included

Unlike the AYMH hardercase, GL Section 6 does not need a numerical internal-mode eigenvalue. It only needs a lower bound excluding eigenvalues below approximately 1.3324. The verifier therefore omits all AYMH blocks for internal-mode Evans functions, threshold Wronskians, outgoing Weyl solutions, and FGR coefficients.
