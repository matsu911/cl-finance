(defpackage cl-finance.model
  (:use common-lisp clf.math)
  (:nicknames clf.model)
  (:import-from clf.math square cumulative-normal-distribution)
  (:export black-scholes))

(defun black-scholes (type S K time sigma r q)
  (let* ((d+ (/ (+ (log (/ S K))
		   (* time (+ r (- q) (/ (square sigma) 2))))
		(* sigma (sqrt time))))
	 (d- (- d+ (* sigma (sqrt time)))))
    (flet ((N (z) (cumulative-normal-distribution z)))
      (flet ((call (S K time r q d+ d-)
	       (- (* S (exp (- (* q time))) (N d+))
		  (* K (exp (- (* r time))) (N d-))))
	     (put (S K time r q d+ d-)
	       (- (* K (exp (- (* r time))) (N (- d-)))
		  (* S (exp (- (* q time))) (N (- d+))))))
	(cond ((or (eq type :CALL)
		   (eq type :C))
	       (call S K time r q d+ d-))
	      ((or (eq type :PUT)
		   (eq type :P))
	       (put S K time r q d+ d-))
	      (t (error "type must be :CALL, :C, :PUT or :P")))))))