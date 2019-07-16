#pragma once

#ifndef clip1_
#define clip1_

// �������������



#include "basic.h"
#include <math.h>
#include <iostream>
using namespace std;

/*bool IsInline(double x, double y, double x1, double y1, double x2, double y2)
{
	double f, dite, eps = 1.0e-5, d1, d2, lmd, k1, k2;
	if ((fabs(x - x1) < 1.0e-10&&fabs(y - y1) < 1.0e-10 || fabs(x - x2) < 1.0e-10&&fabs(y - y2) < 1.0e-10))
		return true;
	else
	{
		dite = x2 - x1;
		if (fabs(dite)<1.e-10)
		{
			k1 = x1;
			f = x - x1;

		}
		else
		{
			k1 = (y2 - y1) / dite;
			f = y - y1 - k1*(x - x1);
			f = f / sqrt(k1*k1 + 1.0);
		}
		if (fabs(x - x1) < 1.0e-10)
			k2 = x1;
		else
			k2 = (y - y1) / (x - x1);

		if (fabs(f) < 3.0e-5 && (k1 - k2) < eps)
		{
			d1 = sqrt((x - x1)*(x - x1) + (y - y1)*(y - y1));
			d2 = sqrt((x - x2)*(x - x2) + (y - y2)*(y - y2));
			lmd = d1 / d2;

			if (fabs(1.0 - x*(1.0 + lmd) / (x1 + lmd*x2)) < 1.0e-6 && fabs(1.0 - y*(1.0 + lmd) / (y1 + lmd*y2)) < 1.0e-6)
				return true;
			else
				return false;
		}
		else
			return false;
	}
}

//���Ƿ���������
int  IsInRegion(double x, double y, Vertex line)
{
	double x1, y1, x2, y2, xpt, ypt;
	double minx, maxx, miny, maxy, eps = 1e-10;
	int n = line.ds;
	int i;
	if (fabs(line.x[n - 1] - line.x[0]) < 1.0e-10 && fabs(line.y[n - 1] - line.y[0]) < 1.0e-10)
		n--;
	minx = maxx = line.x[0];
	miny = maxy = line.y[0];
	for (int i = 1; i < n; i++)
	{
		minx = minx < line.x[i] ? minx : line.x[i];
		maxx = maxx > line.x[i] ? maxx : line.x[i];
		miny = miny < line.y[i] ? miny : line.y[i];
		maxy = maxy > line.y[i] ? maxy : line.y[i];
	}
	if (!(x >= minx - eps&&x <= maxx + eps&&y >= miny - eps&&y <= maxy + eps))
		return -1;

	xpt = maxx + 0.5*(maxx - minx);
	ypt = y;

	double xx1[1], yy1[1];
	int ds = 0, bz = 0;
	for (int i = 0; i < n; i++)
	{
		x1 = line.x[i];
		y1 = line.y[i];
		x2 = line.x[(i + 1) % n];
		y2 = line.y[(i + 1) % n];
		if (IsInline(x, y, x1, y1, x2, y2))
			return 1;
		else
		{
			if (fabs(y1 - y) < eps)
				y1 += 0.5;
			if (fabs(y2 - y) < eps)
				y2 += 0.05;
			bz = CalJoint(x, y, xpt, ypt, x1, y1, x2, y2, xx1, yy1);
			if (bz == 1)
				ds++;
		}

	}
	if (ds % 2 == 1)
		return 2;
	else
		return -1;

}*/

/*�����߶ν���
-1~9�������-1���޽��㣻1���ཻ;2:�˵�A��CD֮�䣻3:�˵�C��AB֮�䣻4:�˵�B��CD֮�䣻5:�˵�D��AB֮�䣻
6��A��C�غϣ�7��A��D�غϣ�8��B��C�غϣ�9��B��D�غϣ�*/
/*int CalJoint(double xa, double ya, double xb, double yb, double xc, double yc, double xd, double yd, double *xx, double *yy)
{
	double fa, fb, fc, fd, cab, ccd, eps = 1.0e-6, qab, qcd;
	double k1, k2, ka, kb, kc, kd;
	cab = (xb - xa)*ya - (yb - ya)*xa;
	ccd = (xd - xc)*yc - (yd - yc)*xc;
	qab = sqrt((yb - ya)*(yb - ya) + (xb - xa)*(xb - xa));
	qcd = sqrt((yd - yc)*(yd - yc) + (xd - xc)*(xd - xc));
	fa = (yd - yc)*xa - (xd - xc)*ya + ccd;
	fb = (yd - yc)*xb - (xd - xc)*yb + ccd;
	fc = (yb - ya)*xc - (xb - xa)*yc + cab;
	fd = (yb - ya)*xd - (xb - xa)*yd + cab;
	fa = fa / qcd;
	fb = fb / qcd;
	fc = fc / qab;
	fd = fd / qab;

	if (fabs(xa - xb) < 1.0e-10)//ֱ��ABб��
		k1 = xa;
	else
		k1 = (yb - ya) / (xb - xa);

	if (fabs(xc - xd) < 1.0e-10)//ֱ��CDб��
		k2 = xc;
	else
		k2 = (yd - yc) / (xd - xc);

	if (fabs(xa - xc) < 1.0e-10)//ֱ��ACб��
		ka = xa;
	else
		ka = (ya - yc) / (xa - xc);

	if (fabs(xb - xc) < 1.0e-10)//ֱ��BCб��
		kb = xb;
	else
		kb = (yb - yc) / (xb - xc);

	if (fabs(xa - xc) < 1.0e-10)//ֱ��CAб��
		kc = xc;
	else
		kc = (yc - ya) / (xc - xa);

	if (fabs(xd - xc) < 1.0e-10)//ֱ��DAб��
		kd = xd;
	else
		kd = (yd - ya) / (xd - xa);

	if (fabs(ka - k2) < eps && fabs(fa) < 3.0e-5)//ֱ��AC��б����CD��ͬ����A����CD��
		fa = 0.0;

	if (fabs(kb - k2) < eps && fabs(fb) < 3.0e-5)//ֱ��BC��б����CD��ͬ����B����CD��
		fb = 0.0;

	if (fabs(kc - k1) < eps && fabs(fc) < 3.0e-5)//ֱ��CA��б����AB��ͬ����C����AB��
		fc = 0.0;

	if (fabs(kd - k1) < eps && fabs(fd) < 3.0e-5)//ֱ��DA��б����AB��ͬ����D����AB��
		fd = 0.0;

	if (fa*fb > 0.0)//A��B��CD�߶�ͬ��
		return -1;

	else if (fa*fb<0.0)//A��B��CD�߶����
	{
		if (fc*fd > 0.0)//C��D��AB�߶�ͬ��
			return -1;
		else if (fc*fd < 0.0)//C��D��AB�߶���࣬AB��CD�ཻ
		{
			double lamd = fabs(fa) / fabs(fb);
			*xx = (xa + lamd*xb) / (1.0 + lamd);
			*yy = (ya + lamd*yb) / (1.0 + lamd);
			return 1;
		}
		else//fc*fd =0.0 ��c��d��һ����AB�߶���
		{
			if (fabs(fc)<eps)//c��AB��
			{
				*xx = xc;
				*yy = yc;
				return 3;
			}
			else if (fabs(fd)<eps)//d��AB��
			{
				*xx = xd;
				*yy = yd;
				return 5;
			}
		}
	}
	else//fa*fb =0.0 ,A��B����һ���������CD�߶���
	{
		if (fc*fd > 0.0)//C��D��AB�߶�ͬ��
			return -1;
		else if (fc*fd < 0.0)//C��D��AB�߶���࣬A��B��һ��CD�߶���
		{
			if (fabs(fa)<eps)//A����CD��
			{
				*xx = xa;
				*yy = ya;
				return 2;
			}
			else if (fabs(fb)<eps)//B����CD��
			{
				*xx = xb;
				*yy = yb;
				return 4;
			}
		}
		else//fc *fd=0.0��C��D��һ���������AB�߶���
		{
			if (fabs(fa)<eps && (!(fabs(fb)<eps)))//A����CD�߶ζ˵��غ�
			{
				if (fabs(fc)<eps && (!(fabs(fd)<eps)))//A��C�غ�
				{
					*xx = xa;
					*yy = ya;
					return 6;
				}
				else if (fabs(fd)<eps && (!(fabs(fc)<eps)))//A��D�غ�
				{
					*xx = xa;
					*yy = ya;
					return 7;
				}
			}
			else if (fabs(fb) < eps && (!(fabs(fa) < eps)))//B����CD�߶ζ˵��غ�
			{
				if (fabs(fc)<eps && (!(fabs(fd)<eps)))//B��C�غ�
				{
					*xx = xb;
					*yy = yb;
					return 8;
				}
				else if (fabs(fd)<eps && (!(fabs(fc)<eps)))//B��D�غ�
				{
					*xx = xb;
					*yy = yb;
					return 9;
				}
			}
			else//AB��CD�غϣ�ֻ��������Ƿ�����һ�߶�֮��
			{
				double d1, d2, lam, fmx, fmy;
				d1 = sqrt((xa - xc)*(xa - xc) + (ya - yc)*(ya - yc));
				d2 = sqrt((xa - xd)*(xa - xd) + (ya - yd)*(ya - yd));
				if (d1<1.0e-8)//a��c�غ�
				{
					*xx = xa;
					*yy = ya;
					return 6;
				}
				else if (d2<1.0e-8)//a��d�غ�
				{
					*xx = xa;
					*yy = ya;
					return 7;
				}
				else
				{
					lam = d1 / d2;
					fmx = xc + lam * xd;
					fmy = yc + lam * yd;
					if (fabs(1.0 - xa*(1.0 + lam) / fmx)<eps && fabs(1.0 - ya*(1.0 + lam) / fmy)<eps)//A����CD֮��
					{
						*xx = xa;
						*yy = ya;
						return 2;
					}
					else
					{
						d2 = sqrt((xc - xd)*(xc - xd) + (yc - yd)*(yc - yd));
						double d3 = sqrt((xd - xb)*(xd - xb) + (yd - yb)*(yd - yb));
						if (d2<1.0e-8)//B��C�����غ�
						{
							*xx = xb;
							*yy = yb;
							return 8;
						}
						else
						{
							lam = d1 / d2;
							fmx = xa + lam*xb;
							fmy = ya + lam*yb;
							if (fabs(1.0 - xc*(1.0 + lam) / fmx) < eps && fabs(1.0 - yc*(1.0 + lam) / fmy) < eps)//C����AB֮��
							{
								*xx = xc;
								*yy = yc;
								return 3;
							}
						}
					}
				}
			}
		}
	}
	return -1;

}*/
//����������
void ConvertLineSort(Vertex *line)
{
	double *x,*y;
	x=new double[line->ds];
	y=new double[line->ds];
	for(int i=0;i<line->ds;i++)
	{	
		x[i]=line->x[line->ds-1-i];
		y[i]=line->y[line->ds-1-i];
	//	cout<<x[i]<<' '<<y[i]<<endl;
	}
	for(int j=0;j<line->ds;j++)
	{	
		line->x[j]=x[j];
		line->y[j]=y[j];
	//	cout<<line->x[j]<<' '<<line->y[j]<<endl;
	}
	delete []x;
	delete []y;
}


//����η����жϣ�>0Ϊ��ʱ�룬��֮��Ϊ˳ʱ�롣��������ж�
double RegionDirection(Vertex line)
{
	double area=0.0,eps=1.0e-10;
	int n=line.ds;
	if(fabs(line.x[n-1]-line.x[0])<eps && fabs(line.y[n-1]-line.y[0])<eps)//��β����
		n--;
	for(int i=0;i<n;i++)
	{
		area=area+(line.x[i]*line.y[(i+1)%n]-line.x[(i+1)%n]*line.y[i]);
	//	cout<<line.x[i]<<' '<<line.y[i]<<endl;

	}
	return 0.5*area;//>0Ϊ��ʱ�룬��֮��Ϊ˳ʱ��
}

int InterPtToPt(InterPoint *pts,InterPoint *pte,InterPoint *pta,InterPoint *ptb,InterPoint *InterPt)
{
	//����ֱ�߶�pts-pte ��ֱ�߶�pt0-pt1�Ľ���InterPt

	double detx, dety, d, yc, ya1, yb1,eps=1.0e-10;
	detx = pts->x - pte->x;
	dety = pts->y - pte->y;
	if (detx==0)
	{
		//����ֱ�߶�pts-pte ��ֱ�߶�pt0-pt1�Ľ���InterPt

		InterPt->parau = ((pta->x - pts->x)*(pte->y - pts->y) - (pta->y - pts->y)*(pte->x - pts->x)) / ((pte->x - pts->x)*(ptb->y - pta->y) - (pte->y - pts->y)*(ptb->x - pta->x));
		InterPt->parav = ((pta->x - pts->x)*(ptb->y - pta->y) - (pta->y - pts->y)*(ptb->x - pta->x)) / ((pte->x - pts->x)*(ptb->y - pta->y) - (pte->y - pts->y)*(ptb->x - pta->x));

		if ((InterPt->parau >= 0 && InterPt->parau<1) && (InterPt->parav >= 0 && InterPt->parav<1))
		{
			//�н���	
			InterPt->x = (1 - InterPt->parav)*pts->x + InterPt->parav*pte->x;
			InterPt->y = (1 - InterPt->parav)*pts->y + InterPt->parav*pte->y;
			return 1;
		}
		else//û�н���
			return 0;

	}
	else {
		d = -dety / detx;
		yc = pts->x*d + pts->y;
		ya1 = yc;
		yb1 = yc;
		double y01, y11;
		y01 = pta->x*d + pta->y;
		y11 = ptb->x*d + ptb->y;
		if ((y01 > yc&&y11 < yc) || (y01 < yc&&y11 > yc))
		{
			InterPt->x = pta->x + (ptb->x - pta->x)*(yc - y01) / (y11 - y01);
			if (pts->x>pte->x)
			{
				if ((InterPt->x >= pte->x && InterPt->x < pts->x))
				{
					InterPt->y = yc - InterPt->x*d;
					InterPt->parau = ((pta->x - pts->x)*(pte->y - pts->y) - (pta->y - pts->y)*(pte->x - pts->x)) / ((pte->x - pts->x)*(ptb->y - pta->y) - (pte->y - pts->y)*(ptb->x - pta->x));
					InterPt->parav = ((pta->x - pts->x)*(ptb->y - pta->y) - (pta->y - pts->y)*(ptb->x - pta->x)) / ((pte->x - pts->x)*(ptb->y - pta->y) - (pte->y - pts->y)*(ptb->x - pta->x));
					return 1;
				}
				else
				{
					return 0;
				}

			}
			else
			{
				if (InterPt->x >= pts->x && InterPt->x < pte->x)
				{
					InterPt->y = yc - InterPt->x*d;
					InterPt->parau = ((pta->x - pts->x)*(pte->y - pts->y) - (pta->y - pts->y)*(pte->x - pts->x)) / ((pte->x - pts->x)*(ptb->y - pta->y) - (pte->y - pts->y)*(ptb->x - pta->x));
					InterPt->parav = ((pta->x - pts->x)*(ptb->y - pta->y) - (pta->y - pts->y)*(ptb->x - pta->x)) / ((pte->x - pts->x)*(ptb->y - pta->y) - (pte->y - pts->y)*(ptb->x - pta->x));
					return 1;
				}
				else
					return 0;

			}

		}
		else if (yc == y01&&yc == y11)
		{
			y01 = pta->x*d + pta->y + eps;
			y11 = ptb->x*d + ptb->y + eps;
			if ((y01 > yc&&y11 < yc) || (y01 < yc&&y11 > yc))
			{
				InterPt->x = pta->x + (ptb->x - pta->x)*(yc - y01) / (y11 - y01);
				if (pts->x>pte->x)
				{
					if ((InterPt->x >= pte->x && InterPt->x < pts->x))
					{
						InterPt->y = yc - InterPt->x*d;
						InterPt->parau = ((pta->x - pts->x)*(pte->y - pts->y) - (pta->y - pts->y)*(pte->x - pts->x)) / ((pte->x - pts->x)*(ptb->y - pta->y) - (pte->y - pts->y)*(ptb->x - pta->x));
						InterPt->parav = ((pta->x - pts->x)*(ptb->y - pta->y) - (pta->y - pts->y)*(ptb->x - pta->x)) / ((pte->x - pts->x)*(ptb->y - pta->y) - (pte->y - pts->y)*(ptb->x - pta->x));
						return 1;
					}
					else
					{
						return 0;
					}

				}
				else
				{
					if (InterPt->x >= pts->x && InterPt->x < pte->x)
					{
						InterPt->y = yc - InterPt->x*d;
						InterPt->parau = ((pta->x - pts->x)*(pte->y - pts->y) - (pta->y - pts->y)*(pte->x - pts->x)) / ((pte->x - pts->x)*(ptb->y - pta->y) - (pte->y - pts->y)*(ptb->x - pta->x));
						InterPt->parav = ((pta->x - pts->x)*(ptb->y - pta->y) - (pta->y - pts->y)*(ptb->x - pta->x)) / ((pte->x - pts->x)*(ptb->y - pta->y) - (pte->y - pts->y)*(ptb->x - pta->x));
						return 1;
					}
					else
						return 0;

				}
			}
			else
			{
				return 0;
			}

		}//û�н���
		else if (yc == y11)
		{
			y11 = ptb->x*d + ptb->y + eps;
			if ((y01 > yc&&y11 < yc) || (y01 < yc&&y11 > yc))
			{
				InterPt->x = pta->x + (ptb->x - pta->x)*(yc - y01) / (y11 - y01);
				if (pts->x>pte->x)
				{
					if ((InterPt->x >= pte->x && InterPt->x < pts->x))
					{
						InterPt->y = yc - InterPt->x*d;
						InterPt->parau = ((pta->x - pts->x)*(pte->y - pts->y) - (pta->y - pts->y)*(pte->x - pts->x)) / ((pte->x - pts->x)*(ptb->y - pta->y) - (pte->y - pts->y)*(ptb->x - pta->x));
						InterPt->parav = ((pta->x - pts->x)*(ptb->y - pta->y) - (pta->y - pts->y)*(ptb->x - pta->x)) / ((pte->x - pts->x)*(ptb->y - pta->y) - (pte->y - pts->y)*(ptb->x - pta->x));
						return 1;
					}
					else
					{
						return 0;
					}

				}
				else
				{
					if (InterPt->x >= pts->x && InterPt->x < pte->x)
					{
						InterPt->y = yc - InterPt->x*d;
						InterPt->parau = ((pta->x - pts->x)*(pte->y - pts->y) - (pta->y - pts->y)*(pte->x - pts->x)) / ((pte->x - pts->x)*(ptb->y - pta->y) - (pte->y - pts->y)*(ptb->x - pta->x));
						InterPt->parav = ((pta->x - pts->x)*(ptb->y - pta->y) - (pta->y - pts->y)*(ptb->x - pta->x)) / ((pte->x - pts->x)*(ptb->y - pta->y) - (pte->y - pts->y)*(ptb->x - pta->x));
						return 1;
					}
					else
						return 0;

				}
			}
			else
			{
				return 0;
			}

		}//û�н���
		else if (yc == y01)
		{
			y01 = ptb->x*d + ptb->y + eps;
			if ((y01 > yc&&y11 < yc) || (y01 < yc&&y11 > yc))
			{
				InterPt->x = pta->x + (ptb->x - pta->x)*(yc - y01) / (y11 - y01);
				if (pts->x>pte->x)
				{
					if ((InterPt->x >= pte->x && InterPt->x < pts->x))
					{
						InterPt->y = yc - InterPt->x*d;
						return 1;
					}
					else
					{
						return 0;
					}

				}
				else
				{
					if (InterPt->x >= pts->x && InterPt->x < pte->x)
					{
						InterPt->y = yc - InterPt->x*d;
						return 1;
					}
					else
						return 0;

				}
			}
			else
			{
				return 0;
			}

		}//û�н���
		else
		{
			return 0;
		}
	}

}
double VectorXVector(InterPoint *pt0,InterPoint *pt1,InterPoint *pt2)
{
	//���������Ĳ��k
	/*   V=(xB-xA)��(y-yA)-(x-xA)��(yB-yA) 
	(pt1.x-pt0.x)  (pt1.y-pt0.y)   0
	(pt2.x-pt1.x)  (pt2.y-pt1.y)   0
	*/	
	return (pt1->x-pt0->x)*(pt2->y-pt1->y)-(pt2->x-pt1->x)*(pt1->y-pt0->y);
}



void Polyline_CutToWA(Vertex MainReg, Vertex ClipReg, OutList &polyline, int &morePolyNum)
{
	/**********************************************************************************
	MainReg���������Σ�ClipReg:����ü�����δ��ڣ�*polyline����������
	*******************************************************************************/

	CircleList1 l1;
	CircleList2 l2;
	//1.���ñ��ü����������
	if (RegionDirection(MainReg) > 0.0)
		ConvertLineSort(&MainReg);

	//2.���òü����������
	if (RegionDirection(ClipReg) > 0.0)
		ConvertLineSort(&ClipReg);

	//3.�����������ü��������������
	l1.CreateCircleList(MainReg);
	l2.CreateCircleList(ClipReg);


	//4.ѭ���Ӳü�����������ȡ��һ���ߣ�ѭ���ͱ��ü������ÿ�����ж��ཻ�����󽻵㣬�������������������εĶ��������У�
	int jd = 0;//������
	InterPoint *pt_0, *pt_1, *pt_11, *ptA, *ptB, *ptA1;
	InterPoint	*interPt;
	LinkList l;//��¼��ǰ���ཻ��ļ��ϣ���������

	int inOroutflag = 0;//0����࣬1���ڲ�
	ptA1 = l2.head->next2;
	ptA = ptA1;
	do
	{
		ptB = ptA->next2;
		//ȡ����ε�һ�������ж����ڲ໹�����
		pt_0 = l1.head->next1;
		if (VectorXVector(ptB, ptA, pt_0)*(-1) > 0)//���
		{
			inOroutflag = 0;
		}
		else
		{
			inOroutflag = 1;//�ڲ�
		}
		//ѭ���ӱ��ü��������ȡ��һ�����㣬�ж��Ƿ��н��㣬���н��㣬�����
		if (!l.IsEmpty())
			l.RemoveAll();//��յ�ǰ�ü��ߵĽ���
		pt_11 = l1.head->next1->next1;
		pt_1 = pt_11;
		do
		{
			//ȡ�����ж����ڲ໹�����
			if (VectorXVector(ptB, ptA, pt_1)*(-1) > 0)//���
			{
				//�����һ������Ҳ����࣬
				if (inOroutflag == 0)
					pt_0 = pt_1;

				else {//��һ���������ڲ⣬������Ƿ��н���,���в���������붥�����У���ȥ��ǰ���㣬�����µ�inOroutflag=0
					//���㽻��

					interPt = new InterPoint;

					if (InterPtToPt(ptA, ptB, pt_0, pt_1, interPt) == 1)
					{
						//�н��㣬���øý����ǳ����ʶ
						interPt->flag = -1;
				
						//�ڵ�ǰ�ü��߽��㼯�в��뽻��(������)
						l.InsertElemAtEnd(interPt);

						//���ü�����ζ��������м��뽻��
						interPt->next1 = pt_0->next1;
						pt_0->next1 = interPt;
						jd++;

					}
					else
						delete	interPt;
					inOroutflag = 0;//�´ζ��������
					pt_0 = pt_1;
				}
			}
			else//�ڲ�
			{
				if (inOroutflag == 1)//�ڲ�
					pt_0 = pt_1;

				else
				{
					//��һ����������࣬���н��㣬���㽻�㣬��������붥�����У��ٲ��붥��
					interPt = new InterPoint;

					if (InterPtToPt(ptA, ptB, pt_0, pt_1, interPt) == 1)
					{
						//�н��㣬���øý����ǽ����ʶ
						interPt->flag = 1;

						//�ڵ�ǰ�ü��߽��㼯�в��뽻��
						l.InsertElemAtEnd(interPt);

						//���ü�����ζ��������м��뽻��
						interPt->next1 = pt_0->next1;
						pt_0->next1 = interPt;
						jd++;

					}
					else
						delete	interPt;
					inOroutflag = 1;//�����µ�inOroutflag=1���ڲ�
					pt_0 = pt_1;
				}
			}
			pt_1 = pt_1->next1;
		} while (pt_1 != pt_11);

		//�ڵ�ǰ�ü��߽��㼯��������Xֵ�����������
		if (!l.IsEmpty())
		{
			InterPoint *nPt, *nPt2;
			l.sort();
			int s = l.GetLength();//�������
			nPt = l.GetData(1);
			nPt2 = l.GetData(s);

			nPt2->next2 = ptB;
			ptA->next2 = nPt;
			l.head->next2 = NULL;
/*			if (ptA->x < ptB->x)
			{
				InterPoint *nPt, *nPt2;
				l.sort1();
				int s = l.GetLength();//�������
				nPt = l.GetData(1);
				nPt2 = l.GetData(s);

				nPt2->next2 = ptB;
				ptA->next2 = nPt;
				l.head->next2 = NULL;
			}
			else if (ptA->x > ptB->x)
			{
				InterPoint *nPt, *nPt2;
				l.sort2();
				int s = l.GetLength();//�������
				nPt = l.GetData(1);
				nPt2 = l.GetData(s);

				nPt2->next2 = ptB;
				ptA->next2 = nPt;
				l.head->next2 = NULL;
			}
			else if (ptA->x = ptB->x)
			{
				if (ptA->y < ptB->y)
				{
					InterPoint *nPt, *nPt2;
					l.sort3();
					int s = l.GetLength();//�������
					nPt = l.GetData(1);
					nPt2 = l.GetData(s);

					nPt2->next2 = ptB;
					ptA->next2 = nPt;
					l.head->next2 = NULL;
				}
				else
				{
					InterPoint *nPt, *nPt2;
					l.sort4();
					int s = l.GetLength();//�������
					nPt = l.GetData(1);
					nPt2 = l.GetData(s);

					nPt2->next2 = ptB;
					ptA->next2 = nPt;
					l.head->next2 = NULL;
				}
			}*/
		}
			ptA = ptB;

		}while (ptA != ptA1);
/*		l1.traverseNode();
		cout << endl;
		l2.traverseNode();
		cout << jd << endl;*/
		if (jd == 0)//�޽��㣬2�������1.�ں���2.����
		{
			return;
		}
			else
			{
			/*ѭ���ӱ��ü�����εĶ��������в��ҽ���,Ȼ�󣬴ӱ��ü�������Ѽ����㣬�������㣬��Ӳü�������Ѽ����㣬ֱ���������㣬
			��������ǳ�ʼ���ҵ��Ľ��㣬������Ѽ��������γ�һ������Σ����ý���ı�ʶ��Ϊ�ǽ��㣬���´���һ�����㿪ʼ�Ѽ�
			*/
				morePolyNum=0;//�����γɵĶ���ε�������=0�������γ�ԭ��ʼ�����
				InterPoint *PF,*PP,*PO;
				PF=l1.head->next1;
				CircleList1 *l3;
				if (PF->flag!=1)
				 {
					do//�Ҵ�һ����ʼ����
					{
						PF=PF->next1;
						if(PF->flag==1)
							break;
					}while(PF!=l1.head->next1);
				}
				PP=PF;
				do{
					if(PP->used==0)
					{
						PO=PP;
						//����һ���µ�����������������ָ��������ͷָ����뵽ָ������Out�������polygon���У�
						l3 = new CircleList1;
						morePolyNum++;
						polyline.InsertNode(l3->head);
						PO->used=1;//��PO��ָ�Ľ���ڵ���Ϊʹ�ù�
						l3->InsertNode(PO->x,PO->y);
						do {
						//����PO��ָ�Ľ���ڵ㿪ʼ��next1ָ���򣩵���һ������ڵ㣨��ΪN1��֮ǰ�ı��ü�����������еĽڵ���뵽���������������󣬲�ʹPOָ��N1

						 do{
								PO=PO->next1;
								PO->used=1;
								l3->InsertNode(PO->x,PO->y);
							}while(PO->flag!=1&&PO->flag!=-1);


							//����PO��ָ�Ľ���ڵ㿪ʼ��next2ָ���򣩵���һ������ڵ㣨��ΪN2��֮ǰ�ı��ü�����������еĽڵ���뵽���������������󣬲�ʹPOָ��N2
						 do{
								PO=PO->next2;
								PO->used = 1;
								l3->InsertNode(PO->x,PO->y);
							}while(PO->flag!=1&&PO->flag!=-1);

						}while(PO!=PP);
						//cout << endl;
					//l3->traverseNode();//��ӡ��������
						PP = PP->next1;
					}
					else //ʹPPָ�򱻲ü�����������е���һ������ڵ㣨����һ���������һ������ڵ㣩
					{
						do
						{
							PP = PP->next1;
							if(PP->flag==1)
								break;

						} while (PP != PF);
					}

				}while(PP!=PF);
				}


	}


#endif