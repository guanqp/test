
  //! Generic Extended %Kalman Filter (EKF) template base class.

  //! \par Usage
  //! "The %Kalman filter is a set of mathematical equations that provides an
  //! efficient computational (recursive) solution of the least-squares method.
  //! The filter is very powerful in several aspects: it supports estimations
  //! of past, present, and even future states, and it can do so even when the
  //! precise nature of the modeled system is unknown." (quoted from [02])
  //! \n
  //! This version of the %Kalman filter is in fact a Variable-Dimension
  //! Extended %Kalman Filter (VDEKF). It supports optimized algorithms 
  //! (translated from Fortran - see [01]), even in the presence of 
  //! correlated process or measurement noise.
  //! \n
  //! To use this template class, you must first inherit from it and implement
  //! some virtual functions. See the example page for more informations. Note
  //! that you can copy freely an \c EKFilter-derived class freely : this can
  //! be useful if you need to branch your filter based on some condition.
  //!
  //! \par Notation
  //! We prefered the notation of [02] : here it is. Assume a state vector
  //! \f$ x \f$ (to estimate) and a non-linear process 
  //! function \f$ f \f$ (to model) that describes the
  //! evolution of this state through time, that is :
  //! \f[ x_k = f \left( x_{k-1}, u_{k-1}, w_{k-1} \right) \f]
  //! where \f$ u \f$ is the (known) input vector fed to the process and 
  //! \f$ w \f$ is the (unknown) process noise vector due to uncertainty
  //! and process modeling errors. Further suppose that the (known) process
  //! noise covariance matrix is : \f[ Q = E \left( w w^T \right) \f]
  //! Now, let's assume a (known) measurement vector \f$ z \f$, which depends 
  //! on the current state \f$ x \f$ in the form of a non-linear function
  //! \f$ h \f$ (to model) : \f[ z_k = h \left( x_k, v_k \right) \f]
  //! where \f$ v \f$ is the (unknown) measurement noise vector with
  //! a (known) covariance matrix : \f[ R = E \left( v v^T \right) \f]
  //! Suppose that we have an estimate of the previous state 
  //! \f$ \hat{x}_{k-1} \f$, called a corrected state or an
  //! <em>a posteriori</em> state estimate. We can build a predicted state
  //! (also called an <em>a priori</em> state estimate) by using \f$ f \f$ :
  //! \f[ \tilde{x}_k = f \left( \hat{x}_{k-1}, u_{k-1}, 0 \right) \f]
  //! since the input is known and the process noise, unknown. With this
  //! predicted state, we can get a predicted measurement vector by
  //! using \f$ h \f$ : \f[ \tilde{z}_k = h \left( \tilde{x}_k, 0 \right) \f]
  //! since the measurement noise is unknown. To obtain a linear
  //! least-squares formulation, we need to linearize those two systems. 
  //! Here are first-order Taylor series centered on \f$ \tilde{x}_k \f$:
  //! \f[ x_k \approx f \left( \hat{x}_{k-1}, u_{k-1}, 0 \right)
  //! + \frac{\partial f}{\partial x} \left( \hat{x}_{k-1}, u_{k-1}, 0 \right)
  //!   \left( \Delta x \right)
  //! + \frac{\partial f}{\partial u} \left( \hat{x}_{k-1}, u_{k-1}, 0 \right)
  //!   \left( \Delta u \right)
  //! + \frac{\partial f}{\partial w} \left( \hat{x}_{k-1}, u_{k-1}, 0 \right)
  //!   \left( \Delta w \right) \f]
  //! \f[ \phantom{x_k} = \tilde{x}_k + A \left( x_{k-1} - \hat{x}_{k-1} 
  //! \right) + W w_{k-1} \f]
  //! We can do the same for the other system :
  //! \f[ z_k \approx h \left( \tilde{x}_k, 0 \right) 
  //! + \frac{\partial h}{\partial x} \left( \tilde{x}_k, 0 \right)
  //!   \left( \Delta x \right)
  //! + \frac{\partial h}{\partial v} \left( \tilde{x}_k, 0 \right)
  //!   \left( \Delta v \right) \f]
  //! \f[ \phantom{z_k} = \tilde{z}_k + H \left( x_k - \tilde{x}_k \right)
  //! + V v_k \f]
  //! The user of this class must derive from it, and implement all the 
  //! functions corresponding to \a A, \a W, \a Q, f, \a H, \a V, \a R
  //! and h.
  //!
  //! \par References
  //! [01] Bierman, G. J. "Factorization Methods for Discrete Sequential
  //! Estimation", Academic Press, 1977. \n
  //! [02] Welch, G. and Bishop, G. "An Introduction to the %Kalman Filter",
  //! http://www.cs.unc.edu/~welch/kalman/kalmanIntro.html
  //!
  //! \par Template parameters
  //! - \c T : Type of elements contained in matrices and vectors. Usually 
  //!          \c float or \c double.
  //! - \c BEG : Starting index of matrices and vectors. Can be either 0 or 1.
  //! - \c OQ : Optimize calculations on \a Q. This can be turned on if \a Q 
  //!           is diagonal.
  //! - \c OVR : Optimize calculations on \a V and \a R. This can be turned on
  //!            if \a V and \a R are both diagonal matrices.
  //! - \c DGB : Debug flag. If \c true, then bound-checking will be performed,
  //!            and \c OutOfBoundError exceptions can be thrown.
  //!
  //! \par Type requirements for T
  //! - \c T must be <b>default constructible</b>.
  //! - \c T must be <b>constructible from</b> \c double.
  //! - \c T must be \b assignable.
  //! - \c T must be <b>equality comparable</b>.
  //! - \c T must be \b serializable.
  //! - \c T must support <b>basic arithmetic operations</b>.
  //! .
  //! This means that, if \c t1, \c t2 are instances of \c T, 
  //! \c op is an arithmetic operator (+ - * /),
  //! \c is is of type
  //! \c istream and \c os is of type \c ostream, the following
  //! expressions must be valid :
  //! - \code T(); T t1; \endcode Default constructor
  //! - \code T(0.0); T t1(1.0); \endcode Constructor from \c double
  //! - \code T t1 = t2; T t1(t2); T(t1); \endcode Copy constructor
  //! - \code t1 op t2 \endcode Arithmetic operation, convertible to \c T
  //! - \code -t1 \endcode Negation operator, convertible to \c T.
  //!       Same as : \code T(0.0) - t1; \endcode
  //! - \code t1 = t2; \endcode Assignment operator
  //! - \code t1 op= t2; \endcode Arithmetic inplace operation.
  //!       Same as : \code t1 = t1 op t2; \endcode
  //! - \code t1 == t2 \endcode Equality comparison, convertible to \c bool
  //! - \code is >> t1; \endcode \c operator>>()
  //! - \code os << t1; \endcode \c operator<<()
  //! 
  //! Finally, note that \c operator>>() and \c operator<<() must be
  //! compatible. Also, \c operator&() must not have been overloaded.