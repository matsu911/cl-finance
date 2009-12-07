(defpackage cl-finance.math-test
  (:use common-lisp 
	clf.math
	clf.unit-test)
  (:nicknames clf.math-test)
  (:import-from clf.math
		bisection
		brent)
  (:import-from clf.unit-test
		deftest
		check)
  (:export test-bisection
	   test-brent))

(in-package cl-finance.math-test)

(macrolet 
    ((check-1 (f &body forms)
	      `(check ,@(loop for accuracy in forms 
			   collect
			     `(close-to 1.0
					(funcall ,f (lambda (x) (- (* x x) 1))
						 0.1 1.5 :accuracy ,accuracy)
					,accuracy)))))
  (deftest test-bisection ()
    (check-1 #'bisection
	     1.0e-4
	     1.0e-6
	     1.0e-8))
  (deftest test-brent ()
    (check-1 #'brent
	     1.0e-4
	     1.0e-6
	     1.0e-8)))
