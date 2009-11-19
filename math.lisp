(defpackage cl-finance.math
  (:use common-lisp)
  (:nicknames clf.math)
  (:export erf))

(in-package cl-finance.math)

(defun square (x)
  (* x x))

(defun sign (x)
  (if (minusp x) -1.0 1.0))

(defun polynomial-value (x &rest coef)
  (loop 
     with (a . l) = coef
     with s = 0
     for i in (reverse l)
     do (setq s (* x (+ s i)))
     finally (return (+ a s))))

(defun erf (x)
  (let ((tiny  least-positive-double-float)
	(one   1.00000000000000000000e+00) 
	;; c (float)0.84506291151
	(erx   8.45062911510467529297e-01) 
	;; Coefficients for approximation to  erf on [0,0.84375]
	(efx   1.28379167095512586316e-01) 
	;; (efx8  1.02703333676410069053e+00) 
	(pp0   1.28379167095512558561e-01) 
	(pp1  -3.25042107247001499370e-01) 
	(pp2  -2.84817495755985104766e-02) 
	(pp3  -5.77027029648944159157e-03) 
	(pp4  -2.37630166566501626084e-05) 
	(qq1   3.97917223959155352819e-01) 
	(qq2   6.50222499887672944485e-02) 
	(qq3   5.08130628187576562776e-03) 
	(qq4   1.32494738004321644526e-04) 
	(qq5  -3.96022827877536812320e-06) 
	;; Coefficients for approximation to  erf  in [0.84375,1.25]
	(pa0  -2.36211856075265944077e-03) 
	(pa1   4.14856118683748331666e-01) 
	(pa2  -3.72207876035701323847e-01) 
	(pa3   3.18346619901161753674e-01) 
	(pa4  -1.10894694282396677476e-01) 
	(pa5   3.54783043256182359371e-02) 
	(pa6  -2.16637559486879084300e-03) 
	(qa1   1.06420880400844228286e-01) 
	(qa2   5.40397917702171048937e-01) 
	(qa3   7.18286544141962662868e-02) 
	(qa4   1.26171219808761642112e-01) 
	(qa5   1.36370839120290507362e-02) 
	(qa6   1.19844998467991074170e-02) 
	;; Coefficients for approximation to  erfc in [1.25,1/0.35]
	(ra0  -9.86494403484714822705e-03) 
	(ra1  -6.93858572707181764372e-01) 
	(ra2  -1.05586262253232909814e+01) 
	(ra3  -6.23753324503260060396e+01) 
	(ra4  -1.62396669462573470355e+02) 
	(ra5  -1.84605092906711035994e+02) 
	(ra6  -8.12874355063065934246e+01) 
	(ra7  -9.81432934416914548592e+00) 
	(sa1   1.96512716674392571292e+01) 
	(sa2   1.37657754143519042600e+02) 
	(sa3   4.34565877475229228821e+02) 
	(sa4   6.45387271733267880336e+02) 
	(sa5   4.29008140027567833386e+02) 
	(sa6   1.08635005541779435134e+02) 
	(sa7   6.57024977031928170135e+00) 
	(sa8  -6.04244152148580987438e-02) 
	;; Coefficients for approximation to  erfc in [1/.35,28]
	(rb0  -9.86494292470009928597e-03) 
	(rb1  -7.99283237680523006574e-01) 
	(rb2  -1.77579549177547519889e+01) 
	(rb3  -1.60636384855821916062e+02) 
	(rb4  -6.37566443368389627722e+02) 
	(rb5  -1.02509513161107724954e+03) 
	(rb6  -4.83519191608651397019e+02) 
	(sb1   3.03380607434824582924e+01) 
	(sb2   3.25792512996573918826e+02) 
	(sb3   1.53672958608443695994e+03) 
	(sb4   3.19985821950859553908e+03) 
	(sb5   2.55305040643316442583e+03) 
	(sb6   4.74528541206955367215e+02) 
	(sb7  -2.24409524465858183362e+01)
	(ax (abs x)))
    (cond ((< ax 0.84375) 
	   (if (< ax 3.7252902984e-09)
	       (* x (1+ efx))
	       (let* ((z (square x))
		      (r (polynomial-value z pp0 pp1 pp2 pp3 pp4))
		      (s (polynomial-value z one qq1 qq2 qq3 qq4 qq5))
		      (y (/ r s)))
		 (* x (1+ y)))))
	  ((< ax 1.25)
	   (let* ((z (1- ax))
		  (P (polynomial-value z pa0 pa1 pa2 pa3 pa4 pa5 pa6))
		  (Q (polynomial-value z one qa1 qa2 qa3 qa4 qa5 qa6)))
	     (* (sign x)
		(+ erx (/ P Q)))))
	  ((>= ax 6)
	   (* (sign x)
	      (- one tiny)))
	  (t 
	   (let ((ss (/ one (square ax)))
		 R
		 S)
	     (if (< ax 2.85714285714285)
		 (setq R (polynomial-value ss ra0 ra1 ra2 ra3 ra4 ra5 ra6 ra7)
		       S (polynomial-value ss one sa1 sa2 sa3 sa4 sa5 sa6 sa7 sa8))
		 (setq R (polynomial-value ss rb0 rb1 rb2 rb3 rb4 rb5 rb6)
		       S (polynomial-value ss one sb1 sb2 sb3 sb4 sb5 sb6 sb7)))
	     (let ((r (exp (- (/ R S) (square ax) 0.5625))))
	       (* (sign x)
		  (1- (/ r ax)))))))))
