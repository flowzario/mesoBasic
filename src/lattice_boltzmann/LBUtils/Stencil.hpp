# ifndef STENCIL_H
# define STENCIL_H

# include <vector>
# include <string>


struct Stencil {

	int nn;
	std::vector<int> exi,eyi,ezi;
	std::vector<double> ex,ey,ez,wa;
	void setStencil(const std::string);

};


# endif  // STENCIL_H

