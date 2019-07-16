#pragma once

#ifndef basic_
#define basic_

 // 声明、定义语句


#include <stdio.h>
#include<iostream>
using namespace std;

typedef struct Vertex{
public:
	int ds;//顶点数
	double *x,*y;
}Vertex;

class  InterPoint{
public:
	double x,y;
	int flag;
	int used;
	double parau,parav;
	InterPoint *next1,*next2;
public:
	InterPoint(){
		used = 0;
		flag=0;
		parau=0;
		parav=0;

	}
};


class CircleList1
{
	private:
    	int length;
    public:
    	InterPoint * head;
    	CircleList1()
    	{
			head = new InterPoint;
    		head->next1 = head;
    		head->x = 0;
			head->y = 0;
    		length = 0;
		}
    	~CircleList1()
    	{
    		delete(head);
		}
		void InsertNode(double x,double y);//在尾部位置插入结点
    	void CreateCircleList(Vertex Any);  //创建单向循环链表 
    	void traverseNode(); 		   //遍历链表
/*		void deleteNode(int n);        //删除位置为n的结点 
		void insertNode(int n,int data);//在制定位置插入结点
		int getLength();               //得到链表长度
		bool isEmpty();                //判断链表是否为空 
*/
};
void CircleList1::InsertNode(double x,double y)
{

	InterPoint *p,*q;	
	p = new InterPoint;
	p->x=x;	
	p->y=y;
	if (head->next1 == head)
	{
		head->next1 = p;
		p->next1 = p;
	}
	else
	{
		q = head->next1->next1;
		if (q == head->next1)
		{
			p->next1 = head->next1;
			q->next1 = p;
		}
		else
		{
			while (q->next1 != head->next1)
				q = q->next1;
			p->next1 = head->next1;
			q->next1 = p;

		}

	}

}

void CircleList1::CreateCircleList(Vertex Any)
{
	InterPoint *p,*q;
	int i,ds;
	ds = Any.ds;
	for(i=0;i<ds;i++)//建立链表区
	{
			if (i==0)
			{
				p = new InterPoint;
				p->x = Any.x[i%ds];
				p->y = Any.y[i%ds];
				head->next1 = p;
				p->next1 = p;
			}
			else 
			{
				p = new InterPoint;
				p->x = Any.x[i%ds];
				p->y = Any.y[i%ds];
				q = head->next1->next1;
				if (q== head->next1)
				{
					p->next1 = head->next1;
					q->next1 = p;
				}
				else
				{
					while (q->next1 != head->next1)
						q = q->next1;
					p->next1 = head->next1;
					q->next1 = p;

				}

			}
	}
}


void CircleList1::traverseNode()  //遍历链表 
{
	InterPoint *p;
	p = head->next1;
	do
	{
		printf("%lf %lf\n", p->x, p->y);
		p=p->next1;

	} while (p!= head->next1);
}

class CircleList2
{
	private:
    	int length;

    public:
    	CircleList2()
    	{
    		head = new InterPoint;
    		head->next2 = head;
    		head->x = 0;
			head->y = 0;
    		length = 0;
		}
    	~CircleList2()
    	{
    		delete(head);
		}
		InterPoint * head;
		void CreateCircleList(Vertex Any);  //创建单向循环链表 
    	void traverseNode(); 		   //遍历链表
};

void CircleList2::CreateCircleList(Vertex Any)
{
	InterPoint *p, *q;
	int i, ds;
	ds = Any.ds;
	for (i = 0; i<ds; i++)//建立链表区
	{
		if (i == 0)
		{
			p = new InterPoint;
			p->x = Any.x[i%ds];
			p->y = Any.y[i%ds];
			head->next2 = p;
			p->next2 = p;
		}
		else
		{
			p = new InterPoint;
			p->x = Any.x[i%ds];
			p->y = Any.y[i%ds];
			q = head->next2->next2;
			if (q == head->next2)
			{
				p->next2 = head->next2;
				q->next2 = p;
			}
			else
			{
				while (q->next2 != head->next2)
					q = q->next2;
				p->next2 = head->next2;
				q->next2 = p;

			}

		}
	}
}

void CircleList2::traverseNode()  //遍历链表 
{
	InterPoint *p;
	p = head->next2;
	do
	{
		printf("%lf %lf\n", p->x, p->y);
		p = p->next2;

	} while (p != head->next2);
}

class Out{
public:
	InterPoint * polygon;
	Out *next;
};

class OutList
{
	private:

    	int length;
    public:
    	Out * head;
    	OutList()
    	{
    		head = new Out;
    		head->next = head;
    		head->polygon=NULL;
    		length = 0;
		}
    	~OutList()
    	{
    		delete(head);
		}
		void InsertNode(InterPoint * polyline);//在尾部位置插入结点
};

void OutList::InsertNode(InterPoint * polyline)
{
	Out *p,*q;	
	p = new Out;
	p->polygon=polyline;
	q=head->next;
	while(q!=head)
		q=q->next;
	p->next=head;
	q->next=p;
}



class LinkList{
public:
	 LinkList();                    //构建一个单链表;
	 ~LinkList();					//销毁一个单链表;
	int GetLength();                //获取线性表长度
	InterPoint * GetData(int i);//获取第i个节点
	void RemoveAll();				//删除所有元素；
	bool IsEmpty();					//判断链表是否为空
	void InsertElemAtEnd(InterPoint * newNode);//在尾部插入节点
     InterPoint * head;              //头结点
	// InterPoint * find_preNode( InterPoint* pNode);//找到前驱节点
	 //void exchange_node(InterPoint* pNode1, InterPoint* pNode2);//交换节点
	 void sort();//升序排序
	 void sort1();//x升序排序
	 void sort2();//x降序排序
	 void sort3();//y升序排序
	 void sort4();//y降序排序
};



LinkList::LinkList()//初始化单链表
{
	head = new InterPoint;
	head->x=0;
	head->y=0;
	head->next2=NULL;
}

LinkList::~LinkList()//销毁单链表
{
	delete head;
}

bool LinkList::IsEmpty()
{
	if(head->next2==NULL)
		return true;
	else
		return false;
}

void LinkList::RemoveAll()
 {
     InterPoint * p = head->next2;
     InterPoint * ptemp=NULL;
     while (head->next2)                    //在头结点的下一个节点逐个删除节点
     {
         ptemp = p;
         p = p->next2;
         delete ptemp;
		 head->next2=p;
    }
}

void LinkList::InsertElemAtEnd(InterPoint * newNode)
{
	InterPoint * p = head;              //定义指针p指向头结点
	if (head->next2 == NULL)
	{            //当头结点为空时，设置newNode为头结点
	       head->next2 = newNode;
			newNode->next2=NULL;
	}  
	else                          //循环直到最后一个节点，将newNode放置在最后
	   {
	       while (p->next2!= NULL)
	            p = p->next2;
	       p->next2 = newNode;
		   newNode->next2=NULL;
	   }
	
}



void LinkList::sort1()
{
	InterPoint * pri, *mid, *tai, *p;
	p = head->next2;
	while (p != NULL)
	{
		mid = head->next2;
		if (head->x > head->next2->x)
		{
			head->next2 = mid->next2;
			mid->next2 = head;
			head = mid;
		}
		pri = head;
		mid = head->next2;
		tai = mid->next2;
		while (mid->next2 != NULL)
		{
			if (mid->x > tai->x)
			{
				pri->next2 = mid->next2;
				mid->next2 = tai->next2;
				tai->next2 = mid;
			}
			pri = pri->next2;
			mid = pri->next2;
			tai = mid->next2;
		}
		p = p->next2;
	}

}

void LinkList::sort2()
{
	InterPoint * pri, *mid, *tai, *p;
	p = head->next2;
	while (p != NULL)
	{
		mid = head->next2;
		if (head->x < head->next2->x)
		{
			head->next2 = mid->next2;
			mid->next2 = head;
			head = mid;
		}
		pri = head;
		mid = head->next2;
		tai = mid->next2;
		while (mid->next2 != NULL)
		{
			if (mid->x < tai->x)
			{
				pri->next2 = mid->next2;
				mid->next2 = tai->next2;
				tai->next2 = mid;
			}
			pri = pri->next2;
			mid = pri->next2;
			tai = mid->next2;
		}
		p = p->next2;
	}

}

void LinkList::sort3()
{
	InterPoint * pri, *mid, *tai, *p;
	p = head->next2;
	while (p != NULL)
	{
		mid = head->next2;
		if (head->y > head->next2->y)
		{
			head->next2 = mid->next2;
			mid->next2 = head;
			head = mid;
		}
		pri = head;
		mid = head->next2;
		tai = mid->next2;
		while (mid->next2 != NULL)
		{
			if (mid->y > tai->y)
			{
				pri->next2 = mid->next2;
				mid->next2 = tai->next2;
				tai->next2 = mid;
			}
			pri = pri->next2;
			mid = pri->next2;
			tai = mid->next2;
		}
		p = p->next2;
	}

}

void LinkList::sort4()
{
	InterPoint * pri, *mid, *tai, *p;
	p = head->next2;
	while (p != NULL)
	{
		mid = head->next2;
		if (head->y < head->next2->y)
		{
			head->next2 = mid->next2;
			mid->next2 = head;
			head = mid;
		}
		pri = head;
		mid = head->next2;
		tai = mid->next2;
		while (mid->next2 != NULL)
		{
			if (mid->y < tai->y)
			{
				pri->next2 = mid->next2;
				mid->next2 = tai->next2;
				tai->next2 = mid;
			}
			pri = pri->next2;
			mid = pri->next2;
			tai = mid->next2;
		}
		p = p->next2;
	}

}
void LinkList::sort()
{
	InterPoint * pri, *mid, *tai, *p;
	p = head->next2;
	while (p != NULL)
	{
		mid = head->next2;
		if (head->parav > head->next2->parav)
		{
			head->next2 = mid->next2;
			mid->next2 = head;
			head = mid;
		}
		pri = head;
		mid = head->next2;
		tai = mid->next2;
		while (mid->next2 != NULL)
		{
			if (mid->parav > tai->parav)
			{
				pri->next2 = mid->next2;
				mid->next2 = tai->next2;
				tai->next2 = mid;
			}
			pri = pri->next2;
			mid = pri->next2;
			tai = mid->next2;
		}
		p = p->next2;
	}

}
//获取单链表的长度
int LinkList::GetLength()
{
    int count = 0;                  //定义count计数
    InterPoint *p = head->next2;           //定义p指向头结点
    do               //当指针的下一个地址不为空时，count+1
   {
        count++;                  
       p = p->next2;                //p指向p的下一个地址
	} while (p != NULL);
    return count;                   //返回count的数据
}


InterPoint * LinkList::GetData(int i)
{
	InterPoint *temp= head;
	int j = 0;
	while (temp&&j < i - 1) {
		temp = temp->next2;
		j++;
	}
	if (!temp || j > i - 1) 
		return NULL;
	
	else 
		return temp->next2;
	
}

#endif