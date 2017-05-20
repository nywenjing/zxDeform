#include "zxsimviewer.h"
#include "zxmath.h"
#include "zxneohookeanmaterial.h"
#include "zxnonlinearfem_forcemodel_sparse.h"
#include "zxtimestepperbackwardeuler.h"
#include "zxcontactconstraint.h"
#include "zxalmfacetofacecontact.h"
#include "zxcubaturegeneratornnhtp.h"
#include "zxcubaturemodelfem.h"

#include "igl/per_vertex_normals.h"
#include <GL/gl.h>

zxSimViewer::zxSimViewer()
{


    m_world = zxSimWorld::create();

    //ica contact
    //init_femsim0();

    //alm contact
    //init_femsim1();

    //hexmesh
    //init_femsim2();

    //cubature generation
    //init_femsim3();

    //subspace simulation
    init_femsim4();

}

void zxSimViewer::init()
{

}

void check_validility(zxBVHNode::Ptr node)
{
    if(node->get_aabb() == nullptr)
        std::cout<<"error"<<std::endl;

    if(node->is_leaf())
        return;
    check_validility(node->get_left());
    check_validility(node->get_right());

}

void zxSimViewer::init_femsim0()
{

    std::string volfilename = "../zxDeform/data/torus.1";
    zxTetrahedralMesh::Ptr vmesh = zxTetrahedralMesh::create(volfilename);

    //use C3D10
    //vmesh->convertToC3D10();

    zxRenderMesh::Ptr rmesh = zxRenderMesh::create(vmesh);
    zxCollisionMesh::Ptr cmesh  = rmesh;
    zxMaterial::Ptr material = zxNeoHookeanMaterial::create(1e4,0.0,1e3);
    zxNonlinearFEM_ForceModel_Sparse::Ptr forcemodel = zxNonlinearFEM_ForceModel_Sparse::create(vmesh,material);
    zxTimeStepperBackwardEuler::Ptr stepper = zxTimeStepperBackwardEuler::create(forcemodel);
    zxBody::Ptr body = zxBody::create();
    body->set_render_mesh(rmesh);
    body->set_collision_mesh(cmesh);
    body->set_stepper(stepper);
    m_world->add_body(body);

    body->get_bvh_tree()->refit();
    check_validility(body->get_bvh_tree()->get_root());

    std::string planeName = "../zxDeform/data/plane.obj";
    zxRenderMesh::Ptr pmesh = zxRenderMesh::create(planeName);
    zxBody::Ptr plane = zxBody::create();
    plane->set_render_mesh(pmesh);
    plane->set_collision_mesh(pmesh);

    m_world->add_body(plane);

    m_rad = 6;

}

void zxSimViewer::init_femsim1()
{
    std::string cubename = "../zxDeform/data/cube.1";
    std::string platename = "../zxDeform/data/plate.1";

    zxSolidMesh::Ptr cubemesh = zxTetrahedralMesh::create(cubename);
    zxSolidMesh::Ptr platemesh = zxTetrahedralMesh::create(platename);

    zxRenderMesh::Ptr cubesurf = zxRenderMesh::create(cubemesh);
    zxRenderMesh::Ptr platesurf = zxRenderMesh::create(platemesh);
    zxMaterial::Ptr material = zxNeoHookeanMaterial::create(1e4,0.33,1e2);
    zxNonlinearFEM_ForceModel_Sparse::Ptr cubeforcemodel = zxNonlinearFEM_ForceModel_Sparse::create(cubemesh,material);
    zxNonlinearFEM_ForceModel_Sparse::Ptr plateforcemodel = zxNonlinearFEM_ForceModel_Sparse::create(platemesh,material);

    cubesurf->save("../zxDeform/data/cube.obj");
    platesurf->save("../zxDeform/data/plate.obj");
    zxALMSimulator::Ptr simulator = zxALMSimulator::create();
    real eps = 1e7;
    zxALMFaceToFaceContact::Ptr contactInterface = zxALMFaceToFaceContact::create(platesurf,cubesurf,eps);

    simulator->add_forcemodel(cubeforcemodel);
    simulator->add_forcemodel(plateforcemodel);
    simulator->add_contactInterface(contactInterface);

    for(size_t i = 0; i < cubemesh->get_num_nodes(); i++)
    {
        zxNode::Ptr node = cubemesh->get_node(i);
        if(node->r0[2] > 2 - zxEPSILON)
           node->m_bc[0] = node->m_bc[1] = node->m_bc[2] = zxFixed;
    }

    for(size_t i = 0; i < platemesh->get_num_nodes(); i++)
    {
        zxNode::Ptr node = platemesh->get_node(i);
        node->rt[2] += 0.01;

        if(node->r0[2] < -0.2 + zxEPSILON)
        {
            node->m_bc[0] = node->m_bc[1] = zxFixed;
            node->m_bc[2] = zxPrescribed;
            node->rl[2] = node->r0[2] + 0.2;

            zxNodalForce::Ptr nforce = zxNodalForce::create();
            nforce->m_node  = node.get();
            nforce->m_force = vec3d(0.0,0.0,1.0);

            //simulator->add_nodal_force(nforce);
        }

        node->r0[2] += 0.01;
    }


    for(double dz = 0.05; dz < 0.6; dz+= 0.2)
    {
        for(size_t i = 0; i < platemesh->get_num_nodes(); i++)
        {
            zxNode::Ptr node = platemesh->get_node(i);

            if(node->m_bc[2] == zxPrescribed)
                node->rl[2] = node->r0[2] + dz;
        }

        simulator->do_simulate();
    }





    zxBody::Ptr cubebody = zxBody::create();
    zxBody::Ptr platebody = zxBody::create();
    cubebody->set_render_mesh(cubesurf);
    platebody->set_render_mesh(platesurf);

    m_world->add_body(cubebody);
    m_world->add_body(platebody);

    m_rad = 6;
}

void zxSimViewer::init_femsim2()
{
    std::string volfilename = "../zxDeform/data/sphereHex.abq";

    //zxAbqReader reader(volfilename);
    zxHexahedralMesh::Ptr vmesh = zxHexahedralMesh::create(volfilename);


    zxRenderMesh::Ptr rmesh = zxRenderMesh::create(vmesh,2);
    zxCollisionMesh::Ptr cmesh  = rmesh;
    zxMaterial::Ptr material = zxNeoHookeanMaterial::create(1e6,0.0,1e3);
    zxNonlinearFEM_ForceModel_Sparse::Ptr forcemodel = zxNonlinearFEM_ForceModel_Sparse::create(vmesh,material);
    zxTimeStepperBackwardEuler::Ptr stepper = zxTimeStepperBackwardEuler::create(forcemodel);
    zxBody::Ptr body = zxBody::create();
    body->set_render_mesh(rmesh);
    body->set_collision_mesh(cmesh);
    body->set_stepper(stepper);
    m_world->add_body(body);

    body->get_bvh_tree()->refit();
    check_validility(body->get_bvh_tree()->get_root());

    std::string planeName = "../zxDeform/data/plane.obj";
    zxRenderMesh::Ptr pmesh = zxRenderMesh::create(planeName);
    zxBody::Ptr plane = zxBody::create();
    plane->set_render_mesh(pmesh);
    plane->set_collision_mesh(pmesh);

    m_world->add_body(plane);

    m_rad = 6;


}

void zxSimViewer::init_femsim3()
{
    std::string volfilename = "../zxDeform/data/cuboid.abq";
    std::string modeuname = "../zxDeform/data/cuboidL.U";
    std::string lamdaname = "../zxDeform/data/cuboid.lamda";

    zxTetrahedralMesh::Ptr vmesh = zxTetrahedralMesh::create(volfilename);
    zxMaterial::Ptr material = zxNeoHookeanMaterial::create(1e6,0.0,1e3);
    zxNonlinearFEM_ForceModel_Sparse::Ptr forcemodel = zxNonlinearFEM_ForceModel_Sparse::create(vmesh,material);

    Eigen::MatrixXd U;
    Eigen::MatrixXd Lamda;

    zx_read_matrix(modeuname,U);
    zx_read_matrix(lamdaname,Lamda);



    Eigen::VectorXd vecLamda = Eigen::Map<Eigen::VectorXd>(Lamda.data(),Lamda.rows());


    zxCubatureModelFEM::Ptr         cubatureFem = zxCubatureModelFEM::create(forcemodel,U,vecLamda);
    zxCubatureGeneratorNNHTP::Ptr   generatorNNHTP = zxCubatureGeneratorNNHTP::create(cubatureFem);
    generatorNNHTP->generateSamples(100,1.0);
    generatorNNHTP->generateCubatures(30);

    m_rad = 6;

}

void zxSimViewer::init_femsim4()
{
    std::string volfilename = "../zxDeform/data/torus.1";
    std::string modeuname = "../zxDeform/data/torusNL.U";
    std::string cubaturefile = "../zxDeform/data/torus.cubature";
    zxTetrahedralMesh::Ptr vmesh = zxTetrahedralMesh::create(volfilename);

    zxRenderMesh::Ptr rmesh = zxRenderMesh::create(vmesh);
    zxCollisionMesh::Ptr cmesh  = rmesh;
    zxMaterial::Ptr material = zxNeoHookeanMaterial::create(1e4,0.0,1e3);
    zxNonlinearFEM_ForceModel_Sparse::Ptr forcemodel = zxNonlinearFEM_ForceModel_Sparse::create(vmesh,material);

    {
        Eigen::MatrixXd U;
        zx_read_matrix(modeuname,U);

        FILE* fp = fopen(cubaturefile.c_str(),"r");

        int nc = 0;
        int id = 0;
        float w = 0;

        fscanf(fp,"%d\n",&nc);

        std::vector<std::pair<int,real>> cubatures(nc);
        for(int i = 0; i < nc; i++)
        {
            fscanf(fp,"%d %f\n",&id,&w);
            cubatures[i].first = id;
            cubatures[i].second = w;
        }

        forcemodel->set_reduced_model(U,cubatures);
        fclose(fp);


    }
    zxTimeStepperBackwardEuler::Ptr stepper = zxTimeStepperBackwardEuler::create(forcemodel);
    zxBody::Ptr body = zxBody::create();
    body->set_render_mesh(rmesh);
    body->set_collision_mesh(cmesh);
    body->set_stepper(stepper);
    m_world->add_body(body);

    body->get_bvh_tree()->refit();
    check_validility(body->get_bvh_tree()->get_root());

    std::string planeName = "../zxDeform/data/plane.obj";
    zxRenderMesh::Ptr pmesh = zxRenderMesh::create(planeName);
    zxBody::Ptr plane = zxBody::create();
    plane->set_render_mesh(pmesh);
    plane->set_collision_mesh(pmesh);

    m_world->add_body(plane);

    m_rad = 6;

}

void zxSimViewer::render()
{
    //m_rad = 6;
    mat4d mv_mat = get_modelview_matrix();//zx_gl_lookAt(m_cen - m_forward * m_rad * 2,m_cen,m_up);
    mat4d proj_mat = get_perspective_matrix();//zx_gl_perspective(60.0,width() * 1.0/ height(),0.1,300);

    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixd(mv_mat.data());
    glMatrixMode(GL_PROJECTION);
    glLoadMatrixd(proj_mat.data());

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    glShadeModel (GL_SMOOTH);
    glClearColor (1.0, 1.0, 1.0, 1.0);

    for(size_t ib = 0; ib < m_world->get_num_bodies(); ib++)
    {
        zxBody::Ptr body = m_world->get_body(ib);
        zxRenderMesh::Ptr r_mesh = body->get_render_mesh();

        Eigen::MatrixXd V,N;
        Eigen::MatrixXi F;
        r_mesh->get_Matrix_Format(V,F);
        igl::per_vertex_normals(V,F,N);

        glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,r_mesh->m_ambient.data());
        glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,r_mesh->m_diffuse.data());
        glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,r_mesh->m_specular.data());
        glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,r_mesh->m_shiness);

        glBegin(GL_TRIANGLES);
        for(int el = 0; el < F.rows(); el++)
            for(int i = 0; i < 3; i++)
            {
                vec3d x = V.row(F(el,i));
                vec3d n = N.row(F(el,i));

                glNormal3dv(n.data());
                glVertex3dv(x.data());
            }

        glEnd();
    }

}

void zxSimViewer::animate()
{
    if(m_isAnimated)
        m_world->do_time_step();
}
