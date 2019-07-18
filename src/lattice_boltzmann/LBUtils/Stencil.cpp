
# include "Stencil.hpp"


void Stencil::setStencil(const std::string stype)
{
	
	if (stype == "D2Q9")
	{
		// lattice velocities defined as:
		// index:  0  1  2  3  4  5  6  7  8
		// ------ --------------------------
		// x:      0 +1 -1  0  0 +1 -1 -1 +1
		// y:      0  0  0 +1 -1 +1 -1 +1 -1
		//
		// 7 3 5
		// 2 0 1
		// 6 4 8
				
		nn = 9;
		
		ex.push_back(0.0);  ey.push_back(0.0);  // 0
		ex.push_back(1.0);  ey.push_back(0.0);  // 1
		ex.push_back(-1.0); ey.push_back(0.0);  // 2
		ex.push_back(0.0);  ey.push_back(1.0);  // 3
		ex.push_back(0.0);  ey.push_back(-1.0); // 4
		ex.push_back(1.0);  ey.push_back(1.0);  // 5
		ex.push_back(-1.0); ey.push_back(-1.0); // 6
		ex.push_back(-1.0); ey.push_back(1.0);  // 7
		ex.push_back(1.0);  ey.push_back(-1.0); // 8
		
		wa.push_back(4.0/9.0);
		wa.push_back(1.0/9.0);
		wa.push_back(1.0/9.0);
		wa.push_back(1.0/9.0);
		wa.push_back(1.0/9.0);
		wa.push_back(1.0/36.0);
		wa.push_back(1.0/36.0);
		wa.push_back(1.0/36.0);
		wa.push_back(1.0/36.0);	
		
		exi.push_back(0);   eyi.push_back(0);
		exi.push_back(1);   eyi.push_back(0);
		exi.push_back(-1);  eyi.push_back(0);
		exi.push_back(0);   eyi.push_back(1);
		exi.push_back(0);   eyi.push_back(-1);
		exi.push_back(1);   eyi.push_back(1);
		exi.push_back(-1);  eyi.push_back(-1);
		exi.push_back(-1);  eyi.push_back(1);
		exi.push_back(1);   eyi.push_back(-1);	
			
	}
	
	if (stype == "D3Q19")
	{
		
		nn = 19;
		
		ex.push_back(0.0);  ey.push_back(0.0);  ez.push_back(0.0);  // 0
		ex.push_back(1.0);  ey.push_back(0.0);  ez.push_back(0.0);  // 1
		ex.push_back(-1.0); ey.push_back(0.0);  ez.push_back(0.0);  // 2
		ex.push_back(0.0);  ey.push_back(1.0);  ez.push_back(0.0);  // 3
		ex.push_back(0.0);  ey.push_back(-1.0); ez.push_back(0.0);  // 4
		ex.push_back(0.0);  ey.push_back(0.0);  ez.push_back(1.0);  // 5
		ex.push_back(0.0);  ey.push_back(0.0);  ez.push_back(-1.0); // 6
		ex.push_back(1.0);  ey.push_back(1.0);  ez.push_back(0.0);  // 7
		ex.push_back(-1.0); ey.push_back(-1.0); ez.push_back(0.0);  // 8
		ex.push_back(1.0);  ey.push_back(0.0);  ez.push_back(1.0);  // 9
		ex.push_back(-1.0); ey.push_back(0.0);  ez.push_back(-1.0); // 10
		ex.push_back(0.0);  ey.push_back(1.0);  ez.push_back(1.0);  // 11
		ex.push_back(0.0);  ey.push_back(-1.0); ez.push_back(-1.0); // 12
		ex.push_back(1.0);  ey.push_back(-1.0); ez.push_back(0.0);  // 13
		ex.push_back(-1.0); ey.push_back(1.0);  ez.push_back(0.0);  // 14
		ex.push_back(1.0);  ey.push_back(0.0);  ez.push_back(-1.0); // 15
		ex.push_back(-1.0); ey.push_back(0.0);  ez.push_back(1.0);  // 16
		ex.push_back(0.0);  ey.push_back(1.0);  ez.push_back(-1.0); // 17
		ex.push_back(0.0);  ey.push_back(-1.0); ez.push_back(1.0);  // 18
		
		wa.push_back(1.0/3.0);  // 0
		wa.push_back(1.0/18.0); // 1
		wa.push_back(1.0/18.0); // 2
		wa.push_back(1.0/18.0); // 3
		wa.push_back(1.0/18.0); // 4
		wa.push_back(1.0/18.0); // 5
		wa.push_back(1.0/18.0); // 6
		wa.push_back(1.0/36.0); // 7 
		wa.push_back(1.0/36.0); // 8
		wa.push_back(1.0/36.0); // 9 
		wa.push_back(1.0/36.0); // 10 
		wa.push_back(1.0/36.0); // 11
		wa.push_back(1.0/36.0); // 12
		wa.push_back(1.0/36.0); // 13
		wa.push_back(1.0/36.0); // 14
		wa.push_back(1.0/36.0); // 15
		wa.push_back(1.0/36.0); // 16
		wa.push_back(1.0/36.0); // 17
		wa.push_back(1.0/36.0); // 18
		
		exi.push_back(0);  eyi.push_back(0);  ezi.push_back(0);  // 0
		exi.push_back(1);  eyi.push_back(0);  ezi.push_back(0);  // 1
		exi.push_back(-1); eyi.push_back(0);  ezi.push_back(0);  // 2
		exi.push_back(0);  eyi.push_back(1);  ezi.push_back(0);  // 3
		exi.push_back(0);  eyi.push_back(-1); ezi.push_back(0);  // 4
		exi.push_back(0);  eyi.push_back(0);  ezi.push_back(1);  // 5
		exi.push_back(0);  eyi.push_back(0);  ezi.push_back(-1); // 6
		exi.push_back(1);  eyi.push_back(1);  ezi.push_back(0);  // 7
		exi.push_back(-1); eyi.push_back(-1); ezi.push_back(0);  // 8
		exi.push_back(1);  eyi.push_back(0);  ezi.push_back(1);  // 9
		exi.push_back(-1); eyi.push_back(0);  ezi.push_back(-1); // 10
		exi.push_back(0);  eyi.push_back(1);  ezi.push_back(1);  // 11
		exi.push_back(0);  eyi.push_back(-1); ezi.push_back(-1); // 12
		exi.push_back(1);  eyi.push_back(-1); ezi.push_back(0);  // 13
		exi.push_back(-1); eyi.push_back(1);  ezi.push_back(0);  // 14
		exi.push_back(1);  eyi.push_back(0);  ezi.push_back(-1); // 15
		exi.push_back(-1); eyi.push_back(0);  ezi.push_back(1);  // 16
		exi.push_back(0);  eyi.push_back(1);  ezi.push_back(-1); // 17
		exi.push_back(0);  eyi.push_back(-1); ezi.push_back(1);  // 18
		
	}
	
}