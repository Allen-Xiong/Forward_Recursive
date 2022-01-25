#ifndef TREE_H_INCLUDED
#define TREE_H_INCLUDED
#include <list>
#include <queue>
#include <iostream>
using namespace std;

template<class T>
class tree{
public:
    virtual void clear()=0;
    virtual bool isEmpty()const=0;
    virtual T root(T flag)const=0;
    virtual T parent(T x,T flag)const=0;
    virtual T child(T x,int i,T flag)const=0;
    virtual void remove(T x,int i)=0;
    virtual void traverse()const=0;
};
//孩子兄弟链表示法 带父指针
template<class T>
class myTree:public tree<T>
{
private:
    struct treeNode
    {
        treeNode *child,*parent,*sibling;   //pointer to child parent and sibling
        T data;
        treeNode():child(nullptr),parent(nullptr),sibling(nullptr){}
        treeNode(T item,treeNode *c=nullptr,treeNode *p=nullptr,treeNode *s=nullptr)
        {
            data=item;
            child=c;
            parent=p;
            sibling=s;
        }
        ~treeNode(){}
    };
    treeNode *treeRoot;     //tree root node
    unsigned int nodenum = 0;
public:
    class iterator
    {
        treeNode* curnode;
        bool empty()const;
    public:
        iterator(treeNode* n = nullptr) :curnode(n) {};
        iterator(const iterator& other) { curnode = other.curnode; }
        ~iterator() {};
        iterator parent();
        iterator child();
        iterator eldersib();
        iterator youngersib();
        bool operator==(const iterator& other);
        bool operator!=(const iterator& other);
        iterator& operator=(const iterator& other);
        T& operator*();
    };
    iterator root();
public:
    myTree(){treeRoot=nullptr;}
    myTree(T x) { treeRoot = new treeNode(x, nullptr, nullptr, nullptr); nodenum += 1; }
    myTree(const myTree& other);
    ~myTree();   
    unsigned int size()const;        //return the number of the nodes
    bool addNode(T father,T son);   //add son node under father node
    bool addRoot(T r);              //add root
    bool isEmpty()const{return (treeRoot==nullptr);}
    T root(T flag)const
    {
        if (treeRoot==nullptr)
            return flag;
        else
            return treeRoot->data;
    }
    int degree(T x)const;
    int generation(T x)const;
    T child(T x,int i,T flag)const;
    void remove(T x,int i);
    T parent(T x,T flag)const;
    bool getNodevec(vector<T*>& nvec);                     //get node vec sort by default.
    void clear();
    void getChainToRoot(T v,list<T> &ch);//v->root
    void preOrder()const;   //preorder
    void midOrder()const;   //midorder
    void postOrder()const;  //postorder
    void levelOrder()const; //levelorder
    void traverse()const;   //traverse
private:
    treeNode* find(T x,treeNode* start)const;
    void remove(treeNode *&start);
    void clear(treeNode *&p);
    void preOrder(treeNode *p)const;
    void midOrder(treeNode *p)const;
    void postOrder(treeNode *p)const;
    void getNodevec(treeNode* p, vector<T*>& nvec);
};
template<class T>
int myTree<T>::generation(T x)const
{
    int n=-1;
    treeNode *p=find(x,treeRoot);
    while (p!=nullptr)
    {
        p=p->parent;
        n++;
    }
    return n;
}
//从节点v反向追踪到根 链ch中最后一个元素是根
template<class T>
void myTree<T>::getChainToRoot(T v,list<T> &ch)
{
    ch.clear();
    treeNode* p=find(v,treeRoot);
    while (p!=nullptr)
    {
        ch.push_back(p->data);
        p=p->parent;
    }
}
template<class T>
int myTree<T>::degree(T x)const
{
    treeNode* p=find(x,treeRoot);
    if (p==nullptr)
        return -1;
    if (p->child==nullptr)
        return 0;
    else
    {
        p=p->child;
        int n=0;
        while (p!=nullptr)
        {
            p=p->sibling;
            ++n;
        }
        return n;
    }
}
template<class T>
bool myTree<T>::addNode(T father,T son)
{
    //无根情况不能增加节点
    if (treeRoot==nullptr)
        return false;
    treeNode* p=find(son,treeRoot);
    //存在与son相同的节点,不能增加节点
    if (p!=nullptr)
        return false;
    p=find(father,treeRoot);
    p->child=new treeNode(son,nullptr,p,p->child);
    nodenum += 1;
    return true;
}
template<class T>
bool myTree<T>::addRoot(T r)
{
    if (treeRoot==nullptr)
    {
        treeRoot=new treeNode(r);
        nodenum += 1;
        return true;
    }
    else
        return false;
}
template<class T>
void myTree<T>::traverse()const
{
    if (treeRoot==nullptr)
    {
        cout<<"No Node!"<<endl;
        return;
    }
    queue<treeNode*> parentque;
    queue<treeNode*> childque;
    parentque.push(treeRoot);
    treeNode *p=treeRoot;
    treeNode *tmp;
    int i=0;
    int m;
    while(!parentque.empty()||!childque.empty())
    {
        m=parentque.size();
        for (int j=0;j<m;++j)
        {
            p=parentque.front();
            parentque.pop();
            //cout<<p->data<<" son node:";
            if (p->child!=nullptr)
            {
                childque.push(p->child);
                tmp=p->child;
                while (tmp!=nullptr)
                {
                    //cout<<tmp->data<<" ";
                    tmp=tmp->sibling;
                }
            }
            cout<<endl;
        }
        m=childque.size();
        for (int j=0;j<m;++j)
        {
            p=childque.front();
            childque.pop();
            while(p!=nullptr)
            {
                parentque.push(p);
                p=p->sibling;
            }
        }
        i++;
    }
}
//层次遍历
template<class T>
void myTree<T>::levelOrder()const
{
    if (treeRoot==nullptr)
    {
        cout<<"No Node!"<<endl;
        return;
    }
    queue<treeNode*> parentque;
    queue<treeNode*> childque;
    cout<<"The tree is \n";
    parentque.push(treeRoot);
    treeNode *p=treeRoot;
    int i=0;
    int m;
    while(!parentque.empty()||!childque.empty())
    {
        cout << "On " << i << " level:";
        m = parentque.size();
        for (int j=0;j<m;++j)
        {
            p=parentque.front();
            parentque.pop();
            cout<<p->data<<" ";
            if (p->child!=nullptr)
                childque.push(p->child);
        }
        m=childque.size();
        for (int j=0;j<m;++j)
        {
            p=childque.front();
            childque.pop();
            while(p!=nullptr)
            {
                parentque.push(p);
                p=p->sibling;
            }
        }
        cout<<endl;
        i++;
    }
}
template<class T>
void myTree<T>::midOrder(myTree<T>::treeNode *p)const
{
    if (p==nullptr)
        return;
    midOrder(p->child);
    cout<<p->data<<" ";
    midOrder(p->sibling);
}
template<class T>
void myTree<T>::postOrder(myTree<T>::treeNode *p)const
{
    if (p==nullptr)
        return;
    postOrder(p->sibling);
    postOrder(p->child);
    cout<<p->data<<" ";
}

template<class T>
void myTree<T>::getNodevec(treeNode* p, vector<T*>& nvec)
{
    if (p == nullptr)
        return;
    else
        nvec.push_back(&p->data);
    getNodevec(p->sibling, nvec);
    getNodevec(p->child, nvec);
    return;
}

template<class T>
void myTree<T>::preOrder(myTree<T>::treeNode *p)const
{
    if (p==nullptr)
        return;
    cout<<p->data<<" ";
    preOrder(p->child);
    preOrder(p->sibling);
}
template<class T>
void myTree<T>::postOrder()const
{
    cout<<"Post order\n";
    postOrder(treeRoot);
    cout<<endl;
}
template<class T>
void myTree<T>::midOrder()const
{
    cout<<"Mid order\n";
    midOrder(treeRoot);
    cout<<endl;
}
template<class T>
void myTree<T>::preOrder()const
{
    cout<<"Pre order\n";
    preOrder(treeRoot);
    cout<<endl;
}
template<class T>
void myTree<T>::clear(myTree<T>::treeNode *&p)
{
    if (p==nullptr)
        return;
    clear(p->child);
    clear(p->sibling);
    delete p;
    p=nullptr;
    nodenum -= 1;
}
template<class T>
void myTree<T>::clear()
{
    clear(treeRoot);
}
template<class T>
typename myTree<T>::iterator myTree<T>::root()
{
    return iterator(treeRoot);
}
template<class T>
myTree<T>::myTree(const myTree<T>& other)
{
    if (other.treeRoot == nullptr)
        treeRoot = nullptr;
    else
    {
        treeRoot = new treeNode(other.treeRoot->data, nullptr, nullptr, nullptr);
        queue<treeNode*> parentque;   //父队列
        parentque.push(other.treeRoot);
        treeNode* p = nullptr, * c = nullptr;
        while (!parentque.empty())
        {
            p = parentque.front();
            parentque.pop();
            c = p->child;
            while (c != nullptr)
            {
                this->addNode(p->data, c->data);
                parentque.push(c);
                c = c->sibling;
            }
        }
    }
}
template<class T>
myTree<T>::~myTree()
{
    clear(treeRoot);
}
template<class T>
unsigned int myTree<T>::size() const
{
    return nodenum;
}
template<class T>
typename myTree<T>::treeNode* myTree<T>::find(T x,myTree<T>::treeNode* start)const
{
    if (start==nullptr) return nullptr;
    if (start->data==x)
        return start;
    treeNode *p=nullptr;
    if (start->child!=nullptr)
        p=find(x,start->child);
    if ((p==nullptr)&&(start->sibling!=nullptr))
        p=find(x,start->sibling);
    return p;
}
template<class T>
T myTree<T>::child(T x,int i,T flag)const
{
    treeNode* p=find(x,treeRoot);
    if (p==nullptr) return flag;
    treeNode* tmp=p->child;
    while (--i>=0&&tmp->sibling!=nullptr)
    {
        tmp=tmp->sibling;
    }
    if (i>=0)
        return flag;
    else
        return tmp->data;

}
template<class T>
void myTree<T>::remove(myTree<T>::treeNode *&start)
{
    //递归实现删除start下的子树包括start
    if (start==nullptr)
    {
        return;
    }
    if (start->child!=nullptr)
        remove(start->child);
    if (start->sibling!=nullptr)
        remove(start->sibling);
    delete start;
    start=nullptr;
    nodenum -= 1;

}
template<class T>
void myTree<T>::remove(T x,int i)
{
    treeNode *p=find(x,treeRoot);
    treeNode *delp=p->child;
    treeNode *pre;
    for (int j=0;j<i;++j)
    {
        if (delp!=nullptr)
            delp=delp->sibling;
        else
        {
            //std::cerr<<"node "<<x<<" do not have "<<i<<"th subtree"<<endl;
            //exit(-1);
            return;
        }
    }
    if (i==0)
    {
        p->child=delp->sibling;
    }
    else
    {
        pre=p->child;
        for (int j=0;j<i-1;++j)
        {
            pre=pre->sibling;
        }
        pre->sibling=delp->sibling;
    }
    delp->sibling=nullptr;
    remove(delp);
}
template<class T>
T myTree<T>::parent(T x,T flag)const
{
    treeNode* p=find(x,treeRoot);
    if (p==nullptr||p->parent==nullptr)
        return flag;
    else
        return p->parent->data;
}
template<class T>
bool myTree<T>::getNodevec(vector<T*>& nvec)
{
    nvec.reserve(nodenum);
    getNodevec(treeRoot, nvec);
    sort(nvec.begin(), nvec.end(), [](T* p1, T* p2) {return (*p1) < (*p2); });
    return true;
}
#endif // TREE_H_INCLUDED

template<class T>
bool myTree<T>::iterator::empty() const
{
    return (curnode == nullptr);
}

template<class T>
typename myTree<T>::iterator myTree<T>::iterator::parent()
{
    if (curnode == nullptr)
        return iterator();
    else
    {
        return iterator(curnode->parent);
    }
}

template<class T>
typename myTree<T>::iterator myTree<T>::iterator::child()
{
    if (curnode != nullptr)
        return iterator(curnode->child);
    else
        return iterator(nullptr);
}

template<class T>
typename myTree<T>::iterator myTree<T>::iterator::eldersib()
{
    if (curnode == nullptr || curnode->parent == nullptr)
        return iterator();
    else
    {
        auto it = this->parent().child();
        if (it == *this)
            return iterator();
        while (true)
        {
            if (it.youngersib() == *this)
                return iterator(it);
            else
                it = it.youngersib();
        }
    }
}

template<class T>
typename myTree<T>::iterator myTree<T>::iterator::youngersib()
{
    if (curnode == nullptr)
        return iterator();
    else
        return iterator(curnode->sibling);
}

template<class T>
bool myTree<T>::iterator::operator==(const myTree<T>::iterator& other)
{
    return (curnode == other.curnode);
}

template<class T>
bool myTree<T>::iterator::operator!=(const myTree<T>::iterator& other)
{
    return !(*this == other);
}

template<class T>
typename myTree<T>::iterator& myTree<T>::iterator::operator=(const iterator& other)
{
    // TODO: insert return statement here
    curnode = other.curnode;
    return *this;
}

template<class T>
T& myTree<T>::iterator::operator*()
{
    // TODO: insert return statement here
    return curnode->data;
}
