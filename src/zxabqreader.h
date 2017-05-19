#ifndef ZXABQREADER_H
#define ZXABQREADER_H

#include "zxsettings.h"

class zxAbqReader
{
public:
    enum ElementType
    {
        C3D4,
        C3D8,
        C3D10,
        S3,
        S6
    };

    class Node
    {
      public:
        vec3d x;
        int   m_id;
    };

    class Element
    {
    public:
        std::vector<int>    m_node_id;
        ElementType m_type;
    };

public:
    zxAbqReader(std::string filename);

public:
    std::vector<Node>       m_nodes;
    std::vector<Element>    m_elements;
};

#endif // ZXABQREADER_H
