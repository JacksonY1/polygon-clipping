#pragma once

#ifndef basic_
#define basic_

 // �������������


#include <stdio.h>
#include<iostream>
using namespace std;

typedef struct Vertex{
public:
	int ds;//������
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
		void InsertNode(double x,double y);//��β��λ�ò�����
    	void CreateCircleList(Vertex Any);  //��������ѭ������ 
    	void traverseNode(); 		   //��������
/*		void deleteNode(int n);        //ɾ��λ��Ϊn�Ľ�� 
		void insertNode(int n,int data);//���ƶ�λ�ò�����
		int getLength();               //�õ�������
		bool isEmpty();                //�ж������Ƿ�Ϊ�� 
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
	for(i=0;i<ds;i++)//����������
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


void CircleList1::traverseNode()  //�������� 
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
		void CreateCircleList(Vertex Any);  //��������ѭ������ 
    	void traverseNode(); 		   //��������
};

void CircleList2::CreateCircleList(Vertex Any)
{
	InterPoint *p, *q;
	int i, ds;
	ds = Any.ds;
	for (i = 0; i<ds; i++)//����������
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

void CircleList2::traverseNode()  //�������� 
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
		void InsertNode(InterPoint * polyline);//��β��λ�ò�����
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
	 LinkList();                    //����һ��������;
	 ~LinkList();					//����һ��������;
	int GetLength();                //��ȡ���Ա���
	InterPoint * GetData(int i);//��ȡ��i���ڵ�
	void RemoveAll();				//ɾ������Ԫ�أ�
	bool IsEmpty();					//�ж������Ƿ�Ϊ��
	void InsertElemAtEnd(InterPoint * newNode);//��β������ڵ�
     InterPoint * head;              //ͷ���
	// InterPoint * find_preNode( InterPoint* pNode);//�ҵ�ǰ���ڵ�
	 //void exchange_node(InterPoint* pNode1, InterPoint* pNode2);//�����ڵ�
	 void sort();//��������
	 void sort1();//x��������
	 void sort2();//x��������
	 void sort3();//y��������
	 void sort4();//y��������
};



LinkList::LinkList()//��ʼ��������
{
	head = new InterPoint;
	head->x=0;
	head->y=0;
	head->next2=NULL;
}

LinkList::~LinkList()//���ٵ�����
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
     while (head->next2)                    //��ͷ������һ���ڵ����ɾ���ڵ�
     {
         ptemp = p;
         p = p->next2;
         delete ptemp;
		 head->next2=p;
    }
}

void LinkList::InsertElemAtEnd(InterPoint * newNode)
{
	InterPoint * p = head;              //����ָ��pָ��ͷ���
	if (head->next2 == NULL)
	{            //��ͷ���Ϊ��ʱ������newNodeΪͷ���
	       head->next2 = newNode;
			newNode->next2=NULL;
	}  
	else                          //ѭ��ֱ�����һ���ڵ㣬��newNode���������
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
//��ȡ������ĳ���
int LinkList::GetLength()
{
    int count = 0;                  //����count����
    InterPoint *p = head->next2;           //����pָ��ͷ���
    do               //��ָ�����һ����ַ��Ϊ��ʱ��count+1
   {
        count++;                  
       p = p->next2;                //pָ��p����һ����ַ
	} while (p != NULL);
    return count;                   //����count������
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