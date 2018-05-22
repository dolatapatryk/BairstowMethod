procedure Bairstow (n               : Integer;
                    a               : vector;
                    mit             : Integer;
                    mincorr,zerodet : Extended;
                    var z,w         : cplxvector;
                    var it,st       : Integer);
{---------------------------------------------------------------------------}
{                                                                           }
{  The procedure Bairstow finds approximate values of complex roots of the  }
{  polynomial a[n]x^n + a[n-1]x^(n-1) + ... + a[1]x + a[0] with real        }
{  coefficients a[i] (i=0,1,...,n) by Bairstow's method.                    }
{  Data:                                                                    }
{    n       - the degree of polynomial,                                    }
{    a       - an array containing coefficients of the polynomial (the      }
{              element a[i] should contain the value of coefficient before  }
{              x^i; i=0,1,...,n),                                           }
{    mit     - maximum number of iterations in the method,                  }
{    mincorr - relative accuracy of the roots.                              }
{  Results:                                                                 }
{    z  - an array of records which contains approximate values of roots    }
{         (the value of each root is equal to z[k].re+i*z[k].im, where i    }
{         denotes the imaginary unit and k=1,2,...,n),                      }
{    w  - an array of records which contains the values of polynomial at    }
{         approximate roots (the value of polynomial at each root is equal  }
{         to w[k].re+i*w[k].im, where i denotes the imaginary unit and      }
{         k=1,2,...,n),                                                     }
{    it - number of iterations.                                             }
{  Other parameters:                                                        }
{    zerodet - the value below which the determinant of a system of linear  }
{              equations, occuring in the method, will be set as 0,         }
{    st      - a variable which within the procedure Bairstow is assigned   }
{              the value of:                                                }
{                1, if n<1 or mit<1 or mincorr<=0 or zerodet<=0,            }
{                2, if during the calculations a system of linear           }
{                   equations, which should be solved in the method, is     }
{                   singular or almost singular,                            }
{                3, if the given number of iterations is too small to       }
{                   obtain sufficient accuracies,                           }
{                0, otherwise.                                              }
{              Note: If st<>0, then it and the elements of arrays z and w   }
{                    are not calculated.                                    }
{  Unlocal identifiers:                                                     }
{    vector     - a type identifier of extended array [q0..qn], where q<=0  }
{                 and qn>=n,                                                }
{    cplxvector - a type identifier of the form:                            }
{                   array [q1..qn] of record                                }
{                                     re,im : Extended                      }
{                                   end;                                    }
{                 where q1<=1 and qn>=n.                                    }
{                                                                           }
{---------------------------------------------------------------------------}
var i,k,n1                       : Integer;
    m,p,pq0,pq1,q,q0,q1,q2,q3,q4 : Extended;
    cond,endpq                   : Boolean;
    b                            : vector;
begin
  if (n<1) or (mit<1) or (mincorr<=0) or (zerodet<=0)
    then st:=1
    else begin
           for i:=0 to n do
             b[n-i]:=a[i];
           st:=0;
           i:=1;
           p:=-1;
           q:=-1;
           n1:=n-1;
           cond:=true;
           it:=0;
           repeat
             if n=1
               then begin
                      z[i].re:=-b[1]/b[0];
                      z[i].im:=0;
                      cond:=false
                    end
               else begin
                      if n=2
                        then begin
                               q:=b[0];
                               p:=b[1]/q;
                               q:=b[2]/q;
                               cond:=false
                             end
                        else begin
                               pq0:=1.0e63;
                               pq1:=pq0;
                               endpq:=true;
                               repeat
                                 it:=it+1;
                                 q0:=0;
                                 q1:=0;
                                 q2:=b[0];
                                 q3:=b[1]-p*q2;
                                 for k:=2 to n do
                                   begin
                                     q4:=b[k]-p*q3-q*q2;
                                     q2:=q2-p*q1-q*q0;
                                     q0:=q1;
                                     q1:=q2;
                                     q2:=q3;
                                     q3:=q4
                                   end;
                                 if abs(q2)+abs(q3)<zerodet
                                   then endpq:=false
                                   else begin
                                          m:=q*q0+p*q1;
                                          q4:=q1*q1+m*q0;
                                          if abs(q4)<zerodet
                                            then st:=2
                                            else begin
                                                   q0:=(q1*q2-q0*q3)/q4;
                                                   q1:=(m*q2+q1*q3)/q4;
                                                   q2:=abs(q0);
                                                   q3:=abs(q1);
                                                   if (q2>mincorr)
                                                     or (q3>mincorr)
                                                     or (q2<pq0)
                                                     and (q3<pq1)
                                                     then begin
                                                            p:=p+q0;
                                                            pq0:=q2;
                                                            q:=q+q1;
                                                            pq1:=q3
                                                          end
                                                     else endpq:=false;
                                                   if (it=mit) and endpq
                                                     then st:=3
                                                 end
                                        end
                               until (st<>0) or not endpq
                             end;
                      if st=0
                        then begin
                               m:=-p/2;
                               q0:=m*m-q;
                               q1:=sqrt(abs(q0));
                               if q0<0
                                 then begin
                                        z[i].re:=m;
                                        z[i+1].re:=m;
                                        z[i].im:=q1;
                                        z[i+1].im:=-q1
                                      end
                                 else begin
                                        if m>0
                                          then m:=m+q1
                                          else m:=m-q1;
                                        z[i+1].re:=m;
                                        if abs(m)=0
                                          then z[i].re:=0
                                          else z[i].re:=q/m;
                                        z[i].im:=0;
                                        z[i+1].im:=0
                                      end;
                               if n>2
                                 then begin
                                        i:=i+2;
                                        n:=n-2;
                                        q0:=0;
                                        q1:=b[0];
                                        for k:=1 to n do
                                          begin
                                            q2:=b[k]-p*q1-q*q0;
                                            b[k]:=q2;
                                            q0:=q1;
                                            q1:=q2
                                          end
                                      end
                             end
                    end
           until (st<>0) or not cond;
           if st=0
             then begin
                    n:=n1+1;
                    for i:=1 to n do
                      begin
                        p:=z[i].re;
                        q:=z[i].im;
                        q1:=a[n];
                        if q=0
                          then begin
                                 for k:=n1 downto 0 do
                                   q1:=q1*p+a[k];
                                 q2:=0
                               end
                          else begin
                                 q0:=sqrt(sqr(p)+sqr(q));
                                 q3:=2*arctan(q/(p+q0));
                                 q4:=n*q3;
                                 q2:=q1*sin(q4);
                                 q1:=q1*cos(q4);
                                 for k:=n1 downto 0 do
                                   begin
                                     p:=a[k];
                                     q4:=k*q3;
                                     q1:=q1*q0+p*cos(q4);
                                     q2:=q2*q0+p*sin(q4)
                                   end
                               end;
                        w[i].re:=q1;
                        w[i].im:=q2
                      end
                  end
         end
end;
