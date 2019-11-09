
int orderint = (int) order;
for(int jj = 0; jj < x_length; jj++)
{
    double y0 = y[jj];
    for(int ii = firstknot[jj]; ii < lastknot[jj]; ii++)
    {
        d[ii-1] += y0*bin(&(knots[-1]),ii,orderint,x[jj]);
    }
}

