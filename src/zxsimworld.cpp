#include "zxsimworld.h"
#include "zxcollisionalgorithm.h"
#include "zxcontactconstraint.h"
#include "zxfrictionconstraint.h"
#include "zxmath.h"
#include "zxcglinearsolver.h"
zxSimWorld::zxSimWorld()
{
    m_dt = 0.01;
    m_margin = 0.01;
    m_friction = 0.1;

    m_col_sceme = ALM;
}

void zxSimWorld::set_collision_response_scheme(zxCollisionResponseScheme scheme)
{
    m_col_sceme = scheme;


}

void  zxSimWorld::do_unconstrained_step()
{
    for(size_t ib = 0; ib < m_bodies.size(); ib++)
    {
        zxBody::Ptr body = m_bodies[ib];

        if(body->is_static()) continue;
        body->do_unconstrained_step(m_dt);
    }
}

void zxSimWorld::do_time_step()
{
    for(size_t ib = 0; ib < m_bodies.size(); ib++)
    {
        zxBody::Ptr body = m_bodies[ib];

        if(body->is_static()) continue;

        if(body->get_stepper()->get_force_model()->is_reduced())
            m_col_sceme = ALM;
    }

    m_dynamic_lcps.clear();
    do_proxy();
    do_unconstrained_step();

    if(m_col_sceme == ALM || m_col_sceme == Penalty)
        do_alm();
    else if(m_col_sceme == LCP)
        do_iac();

    size_t iter = 0;
    while(do_ccd() > 0 && iter < 10)
    {
        if(m_col_sceme == ALM || m_col_sceme == Penalty)
            do_alm();
        else if(m_col_sceme == LCP)
            do_iac();

        iter++;
    }

    if(iter == 10)
        std::cout<<"erro: ccd failed"<<std::endl;

    for(size_t ib = 0; ib < m_bodies.size(); ib++)
    {
        zxBody::Ptr body = m_bodies[ib];

        if(body->is_static()) continue;
        body->get_render_mesh()->update_position();
    }


}


size_t  zxSimWorld::do_ccd()
{
    for(size_t ib = 0; ib < m_bodies.size(); ib++)
    {
        zxBody::Ptr body = m_bodies[ib];

        if(body->is_static()) continue;
        body->get_collision_mesh()->update_position();
        body->get_collision_mesh()->update_aabb(true);
        body->get_bvh_tree()->refit();
    }

    zxFaceFace_CCD::Ptr collider = zxFaceFace_CCD::create();
    collider->set_margin(m_margin);
    for(size_t ib0 = 0; ib0 < m_bodies.size(); ib0++)
        for(size_t ib1 = ib0 + 1; ib1 < m_bodies.size(); ib1++)
        {
            zxBody::Ptr body0 = m_bodies[ib0];
            zxBody::Ptr body1 = m_bodies[ib1];

            if(body0->is_static() && body1->is_static())
                continue;

            zxBVHTree::Ptr tree0 = body0->get_bvh_tree();
            zxBVHTree::Ptr tree1 = body1->get_bvh_tree();

            tree0->get_root()->collid(tree1->get_root(),collider);
        }

    std::list<zxContactPoint::Ptr>&  contacts = collider->get_contacts();

    m_dynamic_lcps.clear();
    for(std::list<zxContactPoint::Ptr>::iterator iter = contacts.begin(); iter != contacts.end(); iter++)
    {
        zxContactPoint::Ptr cp = *iter;
        zxContactConstraint::Ptr cc = zxContactConstraint::create(cp.get());
        m_dynamic_lcps.push_back(cc);
    }
    m_debug_ccd_contact.insert(m_debug_ccd_contact.begin(),contacts.begin(),contacts.end());
    m_debug_ccd_contact.resize(contacts.size());

    std::cout<<"# ccd: "<<m_debug_ccd_contact.size()<<std::endl;

    return contacts.size();
}

size_t  zxSimWorld::do_proxy()
{
    for(size_t ib = 0; ib < m_bodies.size(); ib++)
    {
        zxBody::Ptr body = m_bodies[ib];

        if(body->is_static()) continue;
        body->get_collision_mesh()->update_position();
        body->get_collision_mesh()->update_aabb(false);
        body->get_bvh_tree()->refit();
    }

    {
        zxBody::Ptr body = m_bodies[0];
        zxCollisionMesh::Ptr mesh = body->get_collision_mesh();

        int numC = 0;
        for(size_t i = 0; i < mesh->get_num_verts(); i++)
        {
            zxCollisionMesh::Vert::Ptr v = mesh->get_vert(i);
            if(v->x[1] < -0.39)
                numC++;
        }
        std::cout<<"#C: "<<numC<<std::endl;
        if(numC != 0)
        {
            numC = 0;
        }
    }
    zxFaceFace_Proximity::Ptr collider = zxFaceFace_Proximity::create();
    collider->set_margin(m_margin);
    for(size_t ib0 = 0; ib0 < m_bodies.size(); ib0++)
        for(size_t ib1 = ib0 + 1; ib1 < m_bodies.size(); ib1++)
        {
            zxBody::Ptr body0 = m_bodies[ib0];
            zxBody::Ptr body1 = m_bodies[ib1];

            if(body0->is_static() && body1->is_static())
                continue;

            zxBVHTree::Ptr tree0 = body0->get_bvh_tree();
            zxBVHTree::Ptr tree1 = body1->get_bvh_tree();

            tree0->get_root()->collid(tree1->get_root(),collider);
        }

    std::list<zxContactPoint::Ptr>&  contacts = collider->get_contacts();
    std::list<zxFrictionConstraint::Ptr> frictionC;
    for(std::list<zxContactPoint::Ptr>::iterator iter = contacts.begin(); iter != contacts.end(); iter++)
    {
        zxContactPoint::Ptr cp = *iter;
        zxContactConstraint::Ptr cc = zxContactConstraint::create(cp.get());
        m_dynamic_lcps.push_back(cc);

        vec3d dir0,dir1;
        zx_plane_dir(cp->m_normal,dir0,dir1);
        zxFrictionConstraint::Ptr fc0 = zxFrictionConstraint::create(cc,dir0);
        zxFrictionConstraint::Ptr fc1 = zxFrictionConstraint::create(cc,dir0);
        fc0->set_friction(m_friction);
        fc1->set_friction(m_friction);
        frictionC.push_back(fc0);
        frictionC.push_back(fc1);
    }

    m_dynamic_lcps.insert(m_dynamic_lcps.end(),frictionC.begin(),frictionC.end());

    std::cout<<"# proxy: "<<contacts.size()<<std::endl;




    return contacts.size();

}

void zxSimWorld::do_iac()
{
    std::vector<zxLCPConstraint::Ptr> lcps = m_static_lcps;
    lcps.insert(lcps.end(),m_dynamic_lcps.begin(),m_dynamic_lcps.end());

    size_t n_lcps = lcps.size();

    if(n_lcps == 0)
        return;
    bool converge = false;
    int  maxOuterIter = 40;
    int  maxInnerIter = 10;

    std::for_each(m_bodies.begin(),m_bodies.end(),[](zxBody::Ptr body){if(!body->is_static()) body->reset_iac();});

    for(size_t ic = 0; ic < n_lcps; ic++)
    {
        zxLCPConstraint::Ptr lcp = lcps[ic];
        lcp->init();
    }
    //std::cout<<"______________"<<std::endl;
    while(maxOuterIter-- > 0)
    {
        real norm0 = 0;
        for(size_t ic = 0; ic < n_lcps; ic++)
        {
            zxLCPConstraint::Ptr lcp = lcps[ic];
            real residual = lcp->get_residual();
            if(residual > 0)
                norm0 += residual;
        }

        //std::cout<<"iter: "<<maxOuterIter<<"  iac norm: "<<norm0<<std::endl;

        std::for_each(lcps.begin(),lcps.end(),[](zxLCPConstraint::Ptr lcp){lcp->reset_projected_gs();});

        //for(size_t iter = 0; iter < 10; iter++)
        std::for_each(lcps.begin(),lcps.end(),[](zxLCPConstraint::Ptr lcp){lcp->do_projected_gs_solve();});
        std::for_each(lcps.begin(),lcps.end(),[](zxLCPConstraint::Ptr lcp){lcp->zero_lcp_impulse();});
        std::for_each(lcps.begin(),lcps.end(),[](zxLCPConstraint::Ptr lcp){lcp->set_lcp_impulse();});
        std::for_each(m_bodies.begin(),m_bodies.end(),
                      [](zxBody::Ptr body){if(!body->is_static())body->compute_iac_response();});
    }

    for(size_t ib = 0; ib < m_bodies.size(); ib++)
    {
        zxBody::Ptr body = m_bodies[ib];

        if(body->is_static()) continue;
        body->get_stepper()->update_mesh();
    }

}

class ALM_CG_Functor : public zxCGMultiplier
{
    ZX_MAKE_SHARED_MACO_DEFAULT_CREAT(ALM_CG_Functor)
    public:
        Eigen::VectorXd multiply(Eigen::VectorXd& b);
        void            computeResidual(Eigen::VectorXd& b);
        void            updatePosition(Eigen::VectorXd& x);


public:
    std::vector<zxBody::Ptr> m_bodies;
    std::vector<zxLCPConstraint::Ptr>   m_constraints;
};

void zxSimWorld::do_alm()
{
    std::vector<zxLCPConstraint::Ptr> lcps = m_static_lcps;
    lcps.insert(lcps.end(),m_dynamic_lcps.begin(),m_dynamic_lcps.end());

    size_t n_lcps = lcps.size();

    if(n_lcps == 0)
        return;

    std::for_each(m_bodies.begin(),m_bodies.end(),[](zxBody::Ptr body){if(!body->is_static()) body->reset_iac();});

    for(size_t ic = 0; ic < n_lcps; ic++)
    {
        zxLCPConstraint::Ptr lcp = lcps[ic];
        lcp->init();
        lcp->set_penalty(1e7);
    }

    int numDofs = 0;

    for(size_t i = 0; i < m_bodies.size(); i++)
    {
        zxBody::Ptr body = m_bodies[i];
        if(body->is_static())
            continue;
        numDofs += body->get_stepper()->get_dim();
    }

    Eigen::VectorXd x(numDofs),b(numDofs),dx(numDofs);

    x.setZero();
    b.setZero();
    dx.setZero();

    ALM_CG_Functor::Ptr alm_functor = ALM_CG_Functor::create();
    alm_functor->m_bodies = m_bodies;
    alm_functor->m_constraints = lcps;


    zxCGLinearSolver lsolver(alm_functor);
    lsolver.set_maxIter(1000);
    int iter = 0;
    do
    {
        real tresidual = 0;
        for(size_t ic = 0; ic < n_lcps; ic++)
        {
            zxLCPConstraint::Ptr lcp = lcps[ic];
            tresidual += lcp->get_residual();
        }
        //std::cout<<"residual: "<<tresidual<<std::endl;
        alm_functor->computeResidual(b);
        dx = lsolver.solver(b);
        x += dx;
        alm_functor->updatePosition(x);
        if((iter + 1) % 3 == 0 && m_col_sceme == ALM)
            std::for_each(lcps.begin(),lcps.end(),[](zxLCPConstraint::Ptr lcp){lcp->do_alm_augment();});




    }while(iter++ < 10);


}
void ALM_CG_Functor::computeResidual(Eigen::VectorXd& b)
{
    std::for_each(m_constraints.begin(),m_constraints.end(),[](zxLCPConstraint::Ptr lcp){lcp->zero_lcp_impulse();});
    std::for_each(m_constraints.begin(),m_constraints.end(),[](zxLCPConstraint::Ptr lcp){lcp->set_alm_impulse();});

    int bodyDof = 0;
    for(size_t i = 0; i < m_bodies.size(); i++)
    {
        zxBody::Ptr body = m_bodies[i];
        if(body->is_static())
            continue;
        body->get_stepper()->do_project_alm_impulse(b.data() + bodyDof);
        bodyDof += body->get_stepper()->get_dim();
    }
}

void ALM_CG_Functor::updatePosition(Eigen::VectorXd& x)
{
    int bodyDof = 0;
    for(size_t i = 0; i < m_bodies.size(); i++)
    {
        zxBody::Ptr body = m_bodies[i];
        if(body->is_static())
            continue;
        body->get_stepper()->post_update_alm_constraint(x.data() + bodyDof);
        body->get_stepper()->update_mesh();
        bodyDof += body->get_stepper()->get_dim();
    }

}


Eigen::VectorXd ALM_CG_Functor::multiply(Eigen::VectorXd& b)
{
    Eigen::VectorXd Ab(b.size());

    Ab.setZero();
    int bodyDof = 0;
    for(size_t i = 0; i < m_bodies.size(); i++)
    {
        zxBody::Ptr body = m_bodies[i];
        if(body->is_static())
            continue;

        body->get_stepper()->do_cg_alm_multi(Ab.data() + bodyDof,b.data() + bodyDof);
        bodyDof += body->get_stepper()->get_dim();
    }

    for(size_t i = 0; i < m_constraints.size(); i++)
    {
        zxLCPConstraint::Ptr cons = m_constraints[i];
        cons->compute_alm_Ab();
    }

    bodyDof = 0;
    for(size_t i = 0; i < m_bodies.size(); i++)
    {
        zxBody::Ptr body = m_bodies[i];
        if(body->is_static())
            continue;

        body->get_stepper()->do_project_alm_ab(Ab.data() + bodyDof);
        bodyDof += body->get_stepper()->get_dim();
    }

    return Ab;

}
