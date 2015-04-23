#include "LoopShellShape.h"

namespace voom{

  
  //! constructor function with parrametric coords
  void LoopShellShape::update( const VectorXd & Point, const CornerValences & V)
  {
    _coords = Point;

    // need subdivision
    bool needSubdivide = true;
    if ( V(0) == 6 && V(1)==6 && V(2) == 6 ) needSubdivide = false;

    // sub-division matrix
    SubdivisionMatrix  S= SubdivisionMatrix::Zero(12, _nodes);

    _computeSubdivisionMatrix(V, S, needSubdivide);
    //! convert the parametric coords for sub-patch (irregular patches)
    if( needSubdivide ) _convertParaCoords();
		
    // compute shape functions
    _computeFunctions(S, needSubdivide);

    // compute the 1st derivatives of shape functions
    _computeDerivatives(S, needSubdivide);
	
  }


  void LoopShellShape::_computeFunctions( const SubdivisionMatrix & SDMatrix, bool needSubdivide )
  {
    // set the parametric coords
    const double v = _coords(0);
    const double w = _coords(1);
    const double u = 1 - v - w;

    Matrix<double,1,12> boxSplines;
    // computing ...
    // 12 shape functions for the regular patch element
    //
    boxSplines(0)=(u*u*u*u + 2.0*u*u*u*v)/12.0;

    boxSplines(1)=(u*u*u*u + 2.0*u*u*u*w)/12.0;

    boxSplines(2)=(u*u*u*u + 2.0*u*u*u*w + 6.0*u*u*u*v + 6.0*u*u*v*w +
		   12.0*u*u*v*v + 6.0*u*v*v*w + 6.0*u*v*v*v + 2.0*v*v*v*w +
		   v*v*v*v)/12.0;

    boxSplines(3)=(6.0*u*u*u*u + 24.0*u*u*u*w + 24.0*u*u*w*w + 8.0*u*w*w*w + 
		   w*w*w*w + 24.0*u*u*u*v + 60.0*u*u*v*w + 36.0*u*v*w*w + 
		   6.0*v*w*w*w + 24.0*u*u*v*v + 36.0*u*v*v*w + 12.0*v*v*w*w + 
		   8.0*u*v*v*v + 6.0*v*v*v*w + v*v*v*v)/12.0;

    boxSplines(4)=(u*u*u*u + 6.0*u*u*u*w + 12.0*u*u*w*w + 6.0*u*w*w*w + 
		   w*w*w*w + 2.0*u*u*u*v + 6.0*u*u*v*w + 6.0*u*v*w*w + 
		   2.0*v*w*w*w)/12.0;

    boxSplines(5)=(2.0*u*v*v*v + v*v*v*v)/12.0;

    boxSplines(6)=(u*u*u*u + 6.0*u*u*u*w + 12.0*u*u*w*w + 
		   6.0*u*w*w*w + w*w*w*w + 8.0*u*u*u*v + 36.0*u*u*v*w + 
		   36.0*u*v*w*w + 8.0*v*w*w*w + 24.0*u*u*v*v + 60.0*u*v*v*w + 
		   24.0*v*v*w*w + 24.0*u*v*v*v + 24.0*v*v*v*w + 
		   6.0*v*v*v*v)/12.0;

    boxSplines(7)=(u*u*u*u + 8.0*u*u*u*w + 24.0*u*u*w*w + 24.0*u*w*w*w + 
		   6.0*w*w*w*w + 6.0*u*u*u*v + 36.0*u*u*v*w + 60.0*u*v*w*w + 
		   24.0*v*w*w*w + 12.0*u*u*v*v + 36.0*u*v*v*w + 
		   24.0*v*v*w*w + 6.0*u*v*v*v + 8.0*v*v*v*w + v*v*v*v)/12.0;

    boxSplines(8) =(2.0*u*w*w*w + w*w*w*w)/12.0; 
    boxSplines(9)=(2.0*v*v*v*w + v*v*v*v)/12.0;


    boxSplines(10)=(2.0*u*w*w*w + w*w*w*w + 6.0*u*v*w*w + 6.0*v*w*w*w + 
		    6.0*u*v*v*w + 12.0*v*v*w*w + 2.0*u*v*v*v + 
		    6.0*v*v*v*w + v*v*v*v)/12.0;

    boxSplines(11)=(w*w*w*w + 2.0*v*w*w*w)/12.0;
    
    Matrix<double,1,Dynamic> prod(1,_nodes);
    prod = (boxSplines * SDMatrix);
    for (int i=0; i<_nodes; i++) _N[i] = prod[i];
		
  }



  void LoopShellShape::_computeDerivatives( const SubdivisionMatrix & SDMatrix, bool needSubdivide )
  {
    // set the parametric coords
    const double v = _coords(0);
    const double w = _coords(1);
    const double u = 1 - v - w;

    MatrixXd bsDerivatives(2,12); //rows: directions v,w; cols: shape fns

    // 12 * 2 components of the 1st derivatives of the shape function for the
    // regular element
    bsDerivatives(0,0) = (-6.0*v*u*u - 2.0*u*u*u)/12.0;
    bsDerivatives(1,0) = (-6.0*v*u*u - 4.0*u*u*u)/12.0;

    bsDerivatives(0,1) = (-4.0*u*u*u-6.0*u*u*w)/12.0;
    bsDerivatives(1,1) = (-2.0*u*u*u-6.0*u*u*w)/12.0;

    bsDerivatives(0,2) = (-2.0*v*v*v-6.0*v*v*u
			  + 6.0*v*u*u+2.0*u*u*u)/12.0;
    bsDerivatives(1,2) = (-4.0*v*v*v-18.0*v*v*u
			  - 12.0*v*u*u-2.0*u*u*u
			  - 6.0*v*v*w-12.0*v*u*w
			  - 6.0*u*u*w)/12.0;

    bsDerivatives(0,3) = (-4.0*v*v*v-24.0*v*v*u
			  - 24.0*v*u*u-18.0*v*v*w 
			  - 48.0*v*u*w-12.0*u*u*w
			  - 12.0*v*w*w - 12.0*u*w*w
			  - 2.0*w*w*w)/12.0;

    bsDerivatives(1,3) = (-2.0*v*v*v-12.0*v*v*u
			  - 12.0*v*u*u-12.0*v*v*w
			  - 48.0*v*u*w-24.0*u*u*w
			  - 18.0*v*w*w-24.0*u*w*w
			  - 4.0*w*w*w)/12.0;

    bsDerivatives(0,4) = (-6.0*v*u*u-2.0*u*u*u
			  - 12.0*v*u*w-12.0*u*u*w
			  - 6.0*v*w*w-18.0*u*w*w
			  - 4.0*w*w*w)/12.0;

    bsDerivatives(1,4) = (2.0*u*u*u+6.0*u*u*w
			  - 6.0*u*w*w-2.0*w*w*w)/12.0;

    bsDerivatives(0,5) = (2.0*v*v*v+6.0*v*v*u)/12.0;
    bsDerivatives(1,5) = -v*v*v/6.0;

    bsDerivatives(0,6) = (24.0*v*v*u+24.0*v*u*u
			  + 4.0*u*u*u+12.0*v*v*w
			  + 48.0*v*u*w+18.0*u*u*w
			  + 12.0*v*w*w+12.0*u*w*w
			  + 2.0*w*w*w)/12.0;
  
    bsDerivatives(1,6) = (12.0*v*v*u+12.0*v*u*u
			  + 2.0*u*u*u-12.0*v*v*w
			  + 6.0*u*u*w-12.0*v*w*w
			  - 6.0*u*w*w-2.0*w*w*w)/12.0;

    bsDerivatives(0,7) = (-2.0*v*v*v-6.0*v*v*u
			  + 6.0*v*u*u+2.0*u*u*u
			  - 12.0*v*v*w+12.0*u*u*w
			  - 12.0*v*w*w+12.0*u*w*w)/12.0;

    bsDerivatives(1,7) = (2.0*v*v*v+12.0*v*v*u
			  + 18.0*v*u*u+4.0*u*u*u
			  + 12.0*v*v*w+48.0*v*u*w
			  + 24.0*u*u*w+12.0*v*w*w
			  + 24.0*u*w*w)/12.0;

    bsDerivatives(0,8) = -w*w*w/6.0;
    bsDerivatives(1,8) = (6.0*u*w*w+2.0*w*w*w)/12.0;

    bsDerivatives(0,9) = (4.0*v*v*v+6.0*v*v*w)/12.0;
    bsDerivatives(1,9) = v*v*v/6.0;

    bsDerivatives(0,10) = (2.0*v*v*v+6.0*v*v*u
			   + 12.0*v*v*w+12.0*v*u*w
			   + 18.0*v*w*w+6.0*u*w*w
			   + 4.0*w*w*w)/12.0;

    bsDerivatives(1,10)= (4.0*v*v*v+6.0*v*v*u
			  + 18.0*v*v*w+12.0*v*u*w
			  + 12.0*v*w*w+6.0*u*w*w
			  + 2.0*w*w*w)/12.0;

    bsDerivatives(0,11) = w*w*w/6.0;
    bsDerivatives(1,11) = (6.0*v*w*w+4.0*w*w*w)/12.0;
		
    Matrix<double,2,Dynamic> prod(2,_nodes);
    prod = bsDerivatives * SDMatrix;

    double fact = (needSubdivide)? -2.0: 1.0; //need to multiply by -2 if subdivision

    for (int i=0; i<_nodes; i++){
      _DN[i](0) = fact * prod(0,i);
      _DN[i](1) = fact * prod(1,i);
    }   
    return;
		
  }



  void LoopShellShape::_computeSecondDerivatives( const  SubdivisionMatrix & SDMatrix, bool needSubdivide )
  {
    /* second order derivatives of the box spline shape functions */
    /* der(0, *) derivative with respect to vv                     */
    /* der(1, *) derivative with respect to ww                     */
    /* der(2, *) derivative with respect to vw                     */
    const double v = _coords(0);
    const double w = _coords(1);
    const double u = 1 - v - w;
  
    Matrix<double,3,12> bs2ndDerivatives;

    bs2ndDerivatives(0,0) = v*u;
    bs2ndDerivatives(1,0) = v*u+u*u;
    bs2ndDerivatives(2,0) = (12.0*v*u+6.0*u*u)/12.0;

    bs2ndDerivatives(0,1) = u*u+u*w;
    bs2ndDerivatives(1,1) = u*w;
    bs2ndDerivatives(2,1) = (6.0*u*u+12.0*u*w)/12.0;
             
    bs2ndDerivatives(0,2) = -2.0*v*u;
    bs2ndDerivatives(1,2) = v*v+v*u+v*w+u*w;
    bs2ndDerivatives(2,2) = (6.0*v*v-12.0*v*u
			       -6.0*u*u)/12.0;
             
    bs2ndDerivatives(0,3) = v*v-2.0*u*u
      + v*w-2.0*u*w;
    bs2ndDerivatives(1,3) = -2.0*v*u-2.0*u*u
      + v*w+w*w;
    bs2ndDerivatives(2,3) = (6.0*v*v-12.0*u*u
			       + 24.0*v*w+6.0*w*w)/12.0;
             
    bs2ndDerivatives(0,4) = v*u+v*w+u*w+ w*w;
    bs2ndDerivatives(1,4) = -2.0*u*w;
    bs2ndDerivatives(2,4) = (-6.0*u*u-12.0*u*w 
			       + 6.0*w*w)/12.0;
             
    bs2ndDerivatives(0,5) = v*u;
    bs2ndDerivatives(1,5) = 0.0;
    bs2ndDerivatives(2,5) = -v*v/2.0;
             
    bs2ndDerivatives(0,6) = (-24.0*v*v+12.0*u*u-24.0*v*w
			       + 12.0*u*w)/12.0;
    bs2ndDerivatives(1,6) = (-24.0*v*v-24.0*v*u-24.0*v*w
			       - 24.0*u*w)/12.0;
    bs2ndDerivatives(2,6) = (-12.0*v*v+6.0*u*u-24.0*v*w
			       - 12.0*u*w-6.0*w*w)/12.0;
             
    bs2ndDerivatives(0,7) = -2.0*v*u-2.0*v*w-2.0*u*w- 2.0*w*w;
    bs2ndDerivatives(1,7) = v*u+u*u-2.0*v*w - 2.0*w*w;
    bs2ndDerivatives(2,7) = (-6.0*v*v-12.0*v*u+6.0*u*u 
			       - 24.0*v*w-12.0*w*w)/12.0;
             
    bs2ndDerivatives(0,8) = 0.0;
    bs2ndDerivatives(1,8) = u*w;
    bs2ndDerivatives(2,8) = -w*w/2.0; 
             
    bs2ndDerivatives(0,9) = (12.0*v*v+12.0*v*w)/12.0;
    bs2ndDerivatives(1,9) = 0.0;
    bs2ndDerivatives(2,9) = v*v/2.0;
             
    bs2ndDerivatives(0,10)= (12.0*v*u+12.0*v*w+12.0*u*w
			       + 12.0*w*w)/12.0;
    bs2ndDerivatives(1,10)= v*v+v*u+v*w+u*w;
    bs2ndDerivatives(2,10)= (6.0*v*v+12.0*v*u+24.0*v*w 
			       + 12.0*u*w+6.0*w*w)/12.0;
             
    bs2ndDerivatives(0,11)= 0.0;
    bs2ndDerivatives(1,11)= v*w+w*w;
    bs2ndDerivatives(2,11)= w*w/2.0;

    Matrix<double,3,Dynamic> prod(3,_nodes);
    prod = bs2ndDerivatives * SDMatrix;
    double fact = (needSubdivide)? 4.0: 1.0; //need to multiply by 4 if subdivision

    for(int i=0;i<_nodes;i++){
      _DDN[i](0,0) = fact * prod(0,i);
      _DDN[i](1,1) = fact * prod(1,i);
      _DDN[i](0,1) = fact * prod(2,i);
      _DDN[i](1,0) = fact * prod(2,i);
    }
    	
    return;	
  }


  void LoopShellShape::_computeSubdivisionMatrix( const CornerValences & V, 
						 SubdivisionMatrix & S, 
						 bool needSubdivide)
  {
    //
    //  to understand this part, need the formula of Transformation Matrix 
    //  in the document.
    //
    if ( needSubdivide ){
      int N0 = V(0)-2,  N1 = V(1)-2,  N2 = V(2)-2;
  
      ///////////////////////////////////////////////////////////////////////
      //
      //  if the vertices of elements surrounding one element are overlapped,
      //  need new method to construct the Transformation Matrix.
      //
      assert( N0 != 1 || N1 != 1 || N2 != 1);
      //
      ///////////////////////////////////////////////////////////////////////

      //double w0 = ( 0.625 - sqr( ( 0.375 + 0.25*cos(2.0*M_PI/(N0+2)) ) ) )/(N0+2);
      //double w1 = ( 0.625 - sqr( ( 0.375 + 0.25*cos(2.0*M_PI/(N1+2)) ) ) )/(N1+2);
      //double w2 = ( 0.625 - sqr( ( 0.375 + 0.25*cos(2.0*M_PI/(N2+2)) ) ) )/(N2+2);
      //
      // warren way
      double w0 = 0.375/(N0+2);
      double w1 = 0.375/(N1+2);
      double w2 = 0.375/(N2+2);

      double oneMinusNw0 = 1.0 - ( N0 + 2) * w0;
      double oneMinusNw1 = 1.0 - ( N1 + 2) * w1;
      double oneMinusNw2 = 1.0 - ( N2 + 2) * w2;

      if(false){
	std::cout <<"w0 = " << w0 << std::endl
		  <<"w1 = " << w1 << std::endl
		  <<"w2 = " << w2 << std::endl;
	std::cout <<"oneMinusNw0 = " << oneMinusNw0 << std::endl
		  <<"oneMinusNw1 = " << oneMinusNw1 << std::endl
		  <<"oneMinusNw2 = " << oneMinusNw2 << std::endl;
      }
      const int n = static_cast<int>(_nodes);
      assert(n > 0);
 
      double oneEighth = 1.0/8.0;
      double threeEighths = 3.0/8.0;  

      //
      //  first row
      S(0, 0) = oneEighth;
      S(0, 1) = threeEighths;
      S(0, 3) = threeEighths;
      S(0, 4) = oneEighth;
      // 
      // second row
      S(1, 0) = threeEighths;
      S(1, 1) = oneEighth;
      S(1, 3) = threeEighths;
      S(1, n-1) = oneEighth;
      //
      // third row
      S(2, 0) = w1;
      S(2, 1) = oneMinusNw1;
      for (int col=2;col<=3+N1-1;col++)
	S(2,col) = w1;
      //
      // fourth row
      S(3, 0) = threeEighths;  
      S(3, 1) = threeEighths;  
      S(3, 2) = oneEighth;  
      S(3, 3) = oneEighth;  
      //
      // fifth row
      S(4, 0) = oneMinusNw0;
      S(4, 1) = w0;
      S(4, 2) = w0;
      S(4, 3) = w0;
      for (int col=n-N0+1;col<=n-1;col++)
      S(4, col ) = w0;
      //
      // sixth row
      S(5, 1) = threeEighths;
      S(5, 2) = oneEighth;
      S(5, 3+N1-2) = oneEighth;
      S(5, 3+N1-1) = threeEighths;
      //
      // seventh row
      S(6, 0) = oneEighth;
      S(6, 1) = threeEighths;
      S(6, 2) = threeEighths;
      S(6, 3+N1-1) = oneEighth;
      //
      // eighth row
      S(7, 0) = threeEighths;
      S(7, 1) = oneEighth;
      S(7, 2) = threeEighths;
      S(7, 3+N1+N2-2) = oneEighth;
      //
      // ninth row
      S(8, 0) = threeEighths;
      S(8, 2) = oneEighth;
      S(8, 3+N1+N2-2) = threeEighths;
      if ( 3+N1+N2-1 == N1 + N2 + N0  )
	S(8,3) = 1.0/8.0;
      else if ( 3+N1+N2-1 < N1 + N2 + N0  )		
	S(8, 3+N1+N2-1) = oneEighth;
      else
	std::cout << "never thought about this case ..." << std::endl;
      //
      // tenth row
      S(9, 1) = oneEighth;
      S(9, 2) = threeEighths;
      S(9, 3+N1-1) = threeEighths;
      S(9, 3+N1) = oneEighth;
      //
      // eleventh row
      S(10, 0) = w2;
      S(10, 1) = w2;
      S(10, 2) = oneMinusNw2;
      for (int col = 3+N1-1; col <= 3+N1+N2-2; col++)
	S(10, col) = w2;
      //
      // twelfth row
      S(11, 0) = oneEighth;
      S(11, 2) = threeEighths;
      S(11, 3+N1+N2-3) = oneEighth;
      S(11, 3+N1+N2-2) = threeEighths;
      //
    }
    else{
      S(3 ,0 ) = 1.0;
      S(6 ,1 ) = 1.0;
      S(7 ,2 ) = 1.0;
      S(2 ,3 ) = 1.0;
      S(5 ,4 ) = 1.0;
      S(9 ,5 ) = 1.0;
      S(10 ,6) = 1.0;
      S(11 ,7) = 1.0;
      S(8 ,8 ) = 1.0;
      S(4 ,9 ) = 1.0;
      S(1, 10) = 1.0;
      S(0 ,11) = 1.0;						
    }

    const bool Output_Flag = false;
		
    if( Output_Flag ) {
//       std::cout.precision(4);
//       std::cout.setf(std::ios_base::scientific,std::ios_base::floatfield);
      std::cout << "***** Begin Subdivision Matrix *****"<< std::endl;
      //	      << S << std::endl
      for(int i=0; i<S.rows(); i++) {
	std::cout << "    ";
	for(int j=0; j<S.cols(); j++) {
	  std::cout << S(i,j) << "    "; 
	}
	std::cout << "}" << std::endl;
      }
      std::cout << "*****  End Subdivision Matrix  *****"<< std::endl;
    }	
  }



  void LoopShellShape::_convertParaCoords()
  {
    _coords(0) = 1.0 - 2.0 * _coords(0);
    _coords(1) = 1.0 - 2.0 * _coords(1);	
    //! now, the parametric coords are related to the sub-patch
  }


} // namespace voom
