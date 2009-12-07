(defpackage cl-finance.math-test
  (:use common-lisp clf.math)
  (:nicknames clf.math-test)
  (:export ))

(deftest test-bisection ()
  (check
    (close-to 1.0
	      (clf.math:bisection (lambda (x) (- (* x x) 1))
				  0.1 1.5 :accuracy 1.0e-4)
	      1.0e-4)
    (close-to 1.0
	      (clf.math:bisection (lambda (x) (- (* x x) 1))
				  0.1 1.5 :accuracy 1.0e-6)
	      1.0e-6)
    (close-to 1.0
	      (clf.math:bisection (lambda (x) (- (* x x) 1))
				  0.1 1.5 :accuracy 1.0e-8)
	      1.0e-8)))

(deftest test-brent ()
  (check
    (close-to 1.0
	      (clf.math:brent (lambda (x) (- (* x x) 1))
			      0.1 1.5 :accuracy 1.0e-4)
	      1.0e-4)
    (close-to 1.0
	      (clf.math:brent (lambda (x) (- (* x x) 1))
			      0.1 1.5 :accuracy 1.0e-6)
	      1.0e-6)
    (close-to 1.0
	      (clf.math:brent (lambda (x) (- (* x x) 1))
			      0.1 1.5 :accuracy 1.0e-8)
	      1.0e-8)))