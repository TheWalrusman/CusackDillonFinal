double xLower(const double x, const double y)
{
	return 0 * x*y;
}

double xUpper(const double x, const double y)
{
	return 0 * x*y;//
}

double yLower(const double x, const double y)
{
	return (1 - (4 * (x - 0.5)*(x - 0.5))) + 0 * y;
}

double yUpper(const double x, const double y)
{
	return 0 * x*y;
}