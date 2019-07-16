#pragma once

#ifndef clip1_
#define clip1_

// 声明、定义语句



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

//点是否在区域内
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

/*求两线段交点
-1~9种情况：-1：无交点；1：相交;2:端点A在CD之间；3:端点C在AB之间；4:端点B在CD之间；5:端点D在AB之间；
6：A、C重合；7：A、D重合；8：B、C重合；9：B、D重合；*/
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

	if (fabs(xa - xb) < 1.0e-10)//直线AB斜率
		k1 = xa;
	else
		k1 = (yb - ya) / (xb - xa);

	if (fabs(xc - xd) < 1.0e-10)//直线CD斜率
		k2 = xc;
	else
		k2 = (yd - yc) / (xd - xc);

	if (fabs(xa - xc) < 1.0e-10)//直线AC斜率
		ka = xa;
	else
		ka = (ya - yc) / (xa - xc);

	if (fabs(xb - xc) < 1.0e-10)//直线BC斜率
		kb = xb;
	else
		kb = (yb - yc) / (xb - xc);

	if (fabs(xa - xc) < 1.0e-10)//直线CA斜率
		kc = xc;
	else
		kc = (yc - ya) / (xc - xa);

	if (fabs(xd - xc) < 1.0e-10)//直线DA斜率
		kd = xd;
	else
		kd = (yd - ya) / (xd - xa);

	if (fabs(ka - k2) < eps && fabs(fa) < 3.0e-5)//直线AC的斜率与CD相同，则A点在CD上
		fa = 0.0;

	if (fabs(kb - k2) < eps && fabs(fb) < 3.0e-5)//直线BC的斜率与CD相同，则B点在CD上
		fb = 0.0;

	if (fabs(kc - k1) < eps && fabs(fc) < 3.0e-5)//直线CA的斜率与AB相同，则C点在AB上
		fc = 0.0;

	if (fabs(kd - k1) < eps && fabs(fd) < 3.0e-5)//直线DA的斜率与AB相同，则D点在AB上
		fd = 0.0;

	if (fa*fb > 0.0)//A、B在CD线段同侧
		return -1;

	else if (fa*fb<0.0)//A、B在CD线段异侧
	{
		if (fc*fd > 0.0)//C、D在AB线段同侧
			return -1;
		else if (fc*fd < 0.0)//C、D在AB线段异侧，AB与CD相交
		{
			double lamd = fabs(fa) / fabs(fb);
			*xx = (xa + lamd*xb) / (1.0 + lamd);
			*yy = (ya + lamd*yb) / (1.0 + lamd);
			return 1;
		}
		else//fc*fd =0.0 ，c、d有一点在AB线段上
		{
			if (fabs(fc)<eps)//c在AB上
			{
				*xx = xc;
				*yy = yc;
				return 3;
			}
			else if (fabs(fd)<eps)//d在AB上
			{
				*xx = xd;
				*yy = yd;
				return 5;
			}
		}
	}
	else//fa*fb =0.0 ,A、B点有一点或两点在CD线段上
	{
		if (fc*fd > 0.0)//C、D在AB线段同侧
			return -1;
		else if (fc*fd < 0.0)//C、D在AB线段异侧，A、B有一点CD线段上
		{
			if (fabs(fa)<eps)//A点在CD上
			{
				*xx = xa;
				*yy = ya;
				return 2;
			}
			else if (fabs(fb)<eps)//B点在CD上
			{
				*xx = xb;
				*yy = yb;
				return 4;
			}
		}
		else//fc *fd=0.0，C、D有一点或两点在AB线段上
		{
			if (fabs(fa)<eps && (!(fabs(fb)<eps)))//A点与CD线段端点重合
			{
				if (fabs(fc)<eps && (!(fabs(fd)<eps)))//A、C重合
				{
					*xx = xa;
					*yy = ya;
					return 6;
				}
				else if (fabs(fd)<eps && (!(fabs(fc)<eps)))//A、D重合
				{
					*xx = xa;
					*yy = ya;
					return 7;
				}
			}
			else if (fabs(fb) < eps && (!(fabs(fa) < eps)))//B点与CD线段端点重合
			{
				if (fabs(fc)<eps && (!(fabs(fd)<eps)))//B、C重合
				{
					*xx = xb;
					*yy = yb;
					return 8;
				}
				else if (fabs(fd)<eps && (!(fabs(fc)<eps)))//B、D重合
				{
					*xx = xb;
					*yy = yb;
					return 9;
				}
			}
			else//AB和CD重合，只考虑起点是否在另一线段之间
			{
				double d1, d2, lam, fmx, fmy;
				d1 = sqrt((xa - xc)*(xa - xc) + (ya - yc)*(ya - yc));
				d2 = sqrt((xa - xd)*(xa - xd) + (ya - yd)*(ya - yd));
				if (d1<1.0e-8)//a、c重合
				{
					*xx = xa;
					*yy = ya;
					return 6;
				}
				else if (d2<1.0e-8)//a、d重合
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
					if (fabs(1.0 - xa*(1.0 + lam) / fmx)<eps && fabs(1.0 - ya*(1.0 + lam) / fmy)<eps)//A点在CD之间
					{
						*xx = xa;
						*yy = ya;
						return 2;
					}
					else
					{
						d2 = sqrt((xc - xd)*(xc - xd) + (yc - yd)*(yc - yd));
						double d3 = sqrt((xd - xb)*(xd - xb) + (yd - yb)*(yd - yb));
						if (d2<1.0e-8)//B、C两点重合
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
							if (fabs(1.0 - xc*(1.0 + lam) / fmx) < eps && fabs(1.0 - yc*(1.0 + lam) / fmy) < eps)//C点在AB之间
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
//将顶点逆序
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


//多边形方向判断，>0为逆时针，反之，为顺时针。由面积法判断
double RegionDirection(Vertex line)
{
	double area=0.0,eps=1.0e-10;
	int n=line.ds;
	if(fabs(line.x[n-1]-line.x[0])<eps && fabs(line.y[n-1]-line.y[0])<eps)//首尾相连
		n--;
	for(int i=0;i<n;i++)
	{
		area=area+(line.x[i]*line.y[(i+1)%n]-line.x[(i+1)%n]*line.y[i]);
	//	cout<<line.x[i]<<' '<<line.y[i]<<endl;

	}
	return 0.5*area;//>0为逆时针，反之，为顺时针
}

int InterPtToPt(InterPoint *pts,InterPoint *pte,InterPoint *pta,InterPoint *ptb,InterPoint *InterPt)
{
	//计算直线段pts-pte 和直线段pt0-pt1的交点InterPt

	double detx, dety, d, yc, ya1, yb1,eps=1.0e-10;
	detx = pts->x - pte->x;
	dety = pts->y - pte->y;
	if (detx==0)
	{
		//计算直线段pts-pte 和直线段pt0-pt1的交点InterPt

		InterPt->parau = ((pta->x - pts->x)*(pte->y - pts->y) - (pta->y - pts->y)*(pte->x - pts->x)) / ((pte->x - pts->x)*(ptb->y - pta->y) - (pte->y - pts->y)*(ptb->x - pta->x));
		InterPt->parav = ((pta->x - pts->x)*(ptb->y - pta->y) - (pta->y - pts->y)*(ptb->x - pta->x)) / ((pte->x - pts->x)*(ptb->y - pta->y) - (pte->y - pts->y)*(ptb->x - pta->x));

		if ((InterPt->parau >= 0 && InterPt->parau<1) && (InterPt->parav >= 0 && InterPt->parav<1))
		{
			//有交点	
			InterPt->x = (1 - InterPt->parav)*pts->x + InterPt->parav*pte->x;
			InterPt->y = (1 - InterPt->parav)*pts->y + InterPt->parav*pte->y;
			return 1;
		}
		else//没有交点
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

		}//没有交点
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

		}//没有交点
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

		}//没有交点
		else
		{
			return 0;
		}
	}

}
double VectorXVector(InterPoint *pt0,InterPoint *pt1,InterPoint *pt2)
{
	//计算向量的叉乘k
	/*   V=(xB-xA)·(y-yA)-(x-xA)·(yB-yA) 
	(pt1.x-pt0.x)  (pt1.y-pt0.y)   0
	(pt2.x-pt1.x)  (pt2.y-pt1.y)   0
	*/	
	return (pt1->x-pt0->x)*(pt2->y-pt1->y)-(pt2->x-pt1->x)*(pt1->y-pt0->y);
}



void Polyline_CutToWA(Vertex MainReg, Vertex ClipReg, OutList &polyline, int &morePolyNum)
{
	/**********************************************************************************
	MainReg：任意多边形；ClipReg:任意裁剪多边形窗口；*polyline：输出多边形
	*******************************************************************************/

	CircleList1 l1;
	CircleList2 l2;
	//1.设置被裁剪多边形走向
	if (RegionDirection(MainReg) > 0.0)
		ConvertLineSort(&MainReg);

	//2.设置裁剪多边形走向
	if (RegionDirection(ClipReg) > 0.0)
		ConvertLineSort(&ClipReg);

	//3.建立主、被裁剪多边形区域链表
	l1.CreateCircleList(MainReg);
	l2.CreateCircleList(ClipReg);


	//4.循环从裁剪多边形中逐次取出一条边，循环和被裁剪多边形每条边判断相交，并求交点，并将交点加入两个多边形的顶点序列中，
	int jd = 0;//交点数
	InterPoint *pt_0, *pt_1, *pt_11, *ptA, *ptB, *ptA1;
	InterPoint	*interPt;
	LinkList l;//记录当前边相交点的集合，用于排序

	int inOroutflag = 0;//0，外侧，1，内侧
	ptA1 = l2.head->next2;
	ptA = ptA1;
	do
	{
		ptB = ptA->next2;
		//取多边形第一个顶点判断在内侧还是外侧
		pt_0 = l1.head->next1;
		if (VectorXVector(ptB, ptA, pt_0)*(-1) > 0)//外侧
		{
			inOroutflag = 0;
		}
		else
		{
			inOroutflag = 1;//内侧
		}
		//循环从被裁剪多边形中取出一个顶点，判断是否有交点，如有交点，则加入
		if (!l.IsEmpty())
			l.RemoveAll();//清空当前裁剪边的交点
		pt_11 = l1.head->next1->next1;
		pt_1 = pt_11;
		do
		{
			//取顶点判断在内侧还是外侧
			if (VectorXVector(ptB, ptA, pt_1)*(-1) > 0)//外侧
			{
				//如果上一个顶点也在外侧，
				if (inOroutflag == 0)
					pt_0 = pt_1;

				else {//上一个顶点在内测，则计算是否有交点,如有并将交点插入顶点序列，舍去当前顶点，并设新的inOroutflag=0
					//计算交点

					interPt = new InterPoint;

					if (InterPtToPt(ptA, ptB, pt_0, pt_1, interPt) == 1)
					{
						//有交点，设置该交点是出点标识
						interPt->flag = -1;
				
						//在当前裁剪边交点集中插入交点(有问题)
						l.InsertElemAtEnd(interPt);

						//被裁剪多边形顶点序列中加入交点
						interPt->next1 = pt_0->next1;
						pt_0->next1 = interPt;
						jd++;

					}
					else
						delete	interPt;
					inOroutflag = 0;//下次顶点在外侧
					pt_0 = pt_1;
				}
			}
			else//内侧
			{
				if (inOroutflag == 1)//内侧
					pt_0 = pt_1;

				else
				{
					//上一个顶点在外侧，则有交点，计算交点，将交点插入顶点序列，再插入顶点
					interPt = new InterPoint;

					if (InterPtToPt(ptA, ptB, pt_0, pt_1, interPt) == 1)
					{
						//有交点，设置该交点是进点标识
						interPt->flag = 1;

						//在当前裁剪边交点集中插入交点
						l.InsertElemAtEnd(interPt);

						//被裁剪多边形顶点序列中加入交点
						interPt->next1 = pt_0->next1;
						pt_0->next1 = interPt;
						jd++;

					}
					else
						delete	interPt;
					inOroutflag = 1;//设置新的inOroutflag=1在内测
					pt_0 = pt_1;
				}
			}
			pt_1 = pt_1->next1;
		} while (pt_1 != pt_11);

		//在当前裁剪边交点集排序，依据X值升、降序插入
		if (!l.IsEmpty())
		{
			InterPoint *nPt, *nPt2;
			l.sort();
			int s = l.GetLength();//整体插入
			nPt = l.GetData(1);
			nPt2 = l.GetData(s);

			nPt2->next2 = ptB;
			ptA->next2 = nPt;
			l.head->next2 = NULL;
/*			if (ptA->x < ptB->x)
			{
				InterPoint *nPt, *nPt2;
				l.sort1();
				int s = l.GetLength();//整体插入
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
				int s = l.GetLength();//整体插入
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
					int s = l.GetLength();//整体插入
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
					int s = l.GetLength();//整体插入
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
		if (jd == 0)//无交点，2种情况：1.内含；2.相离
		{
			return;
		}
			else
			{
			/*循环从被裁剪多边形的顶点序列中查找进点,然后，从被裁剪多边形搜集顶点，遇到出点，则从裁剪多边形搜集顶点，直到遇到进点，
			如果进点是初始查找到的进点，则此轮搜集结束，形成一个多边形，将该进点的标识设为非进点，重新从下一个进点开始搜集
			*/
				morePolyNum=0;//额外形成的多边形的数量，=0，正在形成原初始多边形
				InterPoint *PF,*PP,*PO;
				PF=l1.head->next1;
				CircleList1 *l3;
				if (PF->flag!=1)
				 {
					do//找打一个起始进点
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
						//建立一个新的输出多边形链表，并将指向该链表的头指针加入到指针链表Out的最后（在polygon域当中）
						l3 = new CircleList1;
						morePolyNum++;
						polyline.InsertNode(l3->head);
						PO->used=1;//将PO所指的交点节点标记为使用过
						l3->InsertNode(PO->x,PO->y);
						do {
						//将从PO所指的交点节点开始（next1指针域）到下一个交点节点（记为N1）之前的被裁剪多边形链表中的节点加入到输出多边形链表的最后，并使PO指向N1

						 do{
								PO=PO->next1;
								PO->used=1;
								l3->InsertNode(PO->x,PO->y);
							}while(PO->flag!=1&&PO->flag!=-1);


							//将从PO所指的交点节点开始（next2指针域）到下一个交点节点（记为N2）之前的被裁剪多边形链表中的节点加入到输出多边形链表的最后，并使PO指向N2
						 do{
								PO=PO->next2;
								PO->used = 1;
								l3->InsertNode(PO->x,PO->y);
							}while(PO->flag!=1&&PO->flag!=-1);

						}while(PO!=PP);
						//cout << endl;
					//l3->traverseNode();//打印结果多边形
						PP = PP->next1;
					}
					else //使PP指向被裁剪多边形链表中的下一个进点节点（即下一个交点的下一个交点节点）
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