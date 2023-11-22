
int orderint = (int) order;
        
for(int jj = 0; jj < x_length; jj++)
{
    double s = 0;
    for(int ii = firstknot[jj]; ii < lastknot[jj]; ii++)
    {
        s += weights[ii-1]*bin(&(knots[-1]),ii,orderint,x[jj]);
    }
    Sx[jj] = s;
}

