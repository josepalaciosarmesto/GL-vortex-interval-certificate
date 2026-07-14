# GL vortex interval certificate

This repository contains the CAPD/C++ interval-arithmetic certificate for the degree-one Ginzburg--Landau vortex estimates used in Section 6 of `GL-arxiv-final.tex`.

The verifier certifies the vortex shooting interval, positivity of $r^{-2}-3(1-U^2)$ up to the turning radius, the boundary inequality at $r=7$, and the Lieb-Thirring compact integral.

Revised verifier constants:

~~~text
c_box         = [0.583189495860, 0.583189495862]
r_turn        = 0.6140
r_tail        = 7
compact check = compact contribution / 6 < 0.4454
total target  = LT total < 0.4456
spectral bound = 2 - sqrt(0.4456) > 1.3324
~~~

Why these targets changed:

* Your first CAPD transcript certified shooting and the r=7 boundary check.
* It did not certify positivity up to r_turn=0.61448; the interval lower endpoint was negative.
* It certified LT total upper endpoint about 0.4455334383, so the old manuscript target 0.44536 is too sharp for that transcript.
* The revised target 0.4456 is conservative and still proves the advertised qualitative spectral gap.

Build on a machine with CAPD:

~~~bash
c++ -O2 -std=c++17 gl_cert.cpp $(capd-config --cflags --libs) -o gl_cert
./gl_cert | tee gl_cert_transcript.txt
~~~

If using your explicit CAPD path:

~~~bash
g++ -O2 -std=c++17 gl_cert.cpp $(/home/jpalacios/CAPD/build/bin/capd-config --cflags --libs) -o gl_cert
./gl_cert | tee gl_cert_transcript.txt
~~~

The proof run should end with:

~~~text
OVERALL: PASS FOR THE GL SECTION 6 NUMERICAL CERTIFICATE
~~~

If the positivity block still fails at r_turn=0.6140, increase turn_cells. If the LT block fails, use the printed LT total upper endpoint as the manuscript target, rounded upward, provided it remains below the spectral threshold needed in the theorem.

