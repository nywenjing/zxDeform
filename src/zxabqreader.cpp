#include "zxabqreader.h"
#include <fstream>
#include <sstream>
zxAbqReader::zxAbqReader(std::string filename)
{
    std::string fname = filename;
    std::ifstream file;
    file.open(fname);



    std::cout<<filename<<" "<<std::endl;

    assert(file.is_open());

    enum ReadFlag
    {
        node,
        element,
        other
    };

    ReadFlag flag = other;

    int numelementNodes = 0;

    std::list<Node> nodes;
    std::list<Element> elements;

    ElementType etype;

    while(!file.eof())
    {
        char fileLine[1024];
        file.getline(fileLine,1024);

        std::string fileLineString = fileLine;

        if(fileLineString == "")
            continue;

        std::transform(fileLineString.begin(),fileLineString.end(),fileLineString.begin(),::tolower);

        if(fileLineString.at(0) == '*')
        {
            flag = other;
            if(fileLineString.find("node") != std::string::npos)
            {
                flag = node;

            }
            else if(fileLineString.find("element") != std::string::npos)
            {
                size_t typeId = fileLineString.find("type=") + 5;
                size_t endTyepId = fileLineString.find(",",typeId);

                std::string type = fileLineString.substr(typeId,endTyepId - typeId);

                if(type == "c3d4")
                {
                    numelementNodes = 4;
                    etype = C3D4;

                }
                else if(type == "c3d10")
                {
                    numelementNodes = 10;
                    etype = C3D10;

                }
                else if(type == "c3d8")
                {
                    numelementNodes = 8;
                    etype = C3D8;
                }

                flag = element;
            }
        }
        else if(flag == node)
        {
            std::replace(fileLineString.begin(),fileLineString.end(),',',' ');
            std::stringstream s(fileLineString);
            int id;
            vec3d x;
            s>>id >>x[0] >> x[1] >>x[2];
            x *= 1.0;

            Node n;
            n.m_id = id;
            n.x = x;
            nodes.push_back(n);
        }
        else if(flag == element)
        {
            std::replace(fileLineString.begin(),fileLineString.end(),',',' ');
            std::stringstream s(fileLineString);
            int e_id;
            int n_id;
            s>>e_id;

            Element el;
            el.m_node_id.resize(numelementNodes);
            for(int i = 0; i < numelementNodes; i++)
            {
                s>>n_id;
                el.m_node_id[i] = n_id;
            }

            elements.push_back(el);
        }

    }

    file.close();

    m_nodes.insert(m_nodes.begin(),nodes.begin(),nodes.end());
    m_elements.insert(m_elements.begin(),elements.begin(),elements.end());
}
