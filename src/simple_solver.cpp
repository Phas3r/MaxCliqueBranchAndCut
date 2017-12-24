#include "simple_solver.h"
#include <algorithm>
#include <iterator>
#include <set>
#include <cassert>

// Overloading operators for set intersections and unions
template <class T, class CMP = std::less<T>, class ALLOC = std::allocator<T> >
std::set<T, CMP, ALLOC> operator * (
    const std::set<T, CMP, ALLOC> &s1, const std::set<T, CMP, ALLOC> &s2)
{
    std::set<T, CMP, ALLOC> s;
    std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
        std::inserter(s, s.begin()));
    return s;
}

template <class T, class CMP = std::less<T>, class ALLOC = std::allocator<T> >
std::set<T, CMP, ALLOC> operator + (
    const std::set<T, CMP, ALLOC> &s1, const std::set<T, CMP, ALLOC> &s2)
{
    std::set<T, CMP, ALLOC> s;
    std::set_union(s1.begin(), s1.end(), s2.begin(), s2.end(),
        std::inserter(s, s.begin()));
    return s;
}


bool validateColoring(const Graph& graph, GraphColoring& coloring)
{
    for (Color c = 0; c < coloring.getColorCount(); ++c) {
        for (Vertex v1 : coloring.getVerticesByColor(c)) {
            for (Vertex v2 : coloring.getVerticesByColor(c)) {
                if (graph.HasEdge(v1, v2)) {
                    std::cout << "UJAS: EDGE IN COLOR" << c << ":" << v1 << v2;
                    assert(false);
                    return false;
                }
            }
        };
    }
    for (Vertex v = 0; v < graph.GetVertexCount(); ++v) {
        if (coloring.getVertexColor(v) >= coloring.getColorCount()) {
            std::cout << "UJAS: UNCOLORED VERTEX" << v;
            assert(false);
            return false;
        };
    }
    return true;
}

void checkClique(Vertices & c, const Graph& g)
{
    for (auto v : c)
    {
        for (auto v2 : c)
        {
            if (v2 != v && !g.HasEdge(v, v2))
                std::cout << "UJAS!"<<std::endl;
        }

    }
}

static bool almostEqual(double x, double y)
{
    double maxXYOne = std::max({ 1.0, std::fabs(x) , std::fabs(y) });

    return std::fabs(x - y) <= std::numeric_limits<double>::epsilon()*maxXYOne;
}

static bool isIntegerVar(double var)
{
    return almostEqual(var, 1.) || almostEqual(var, 0.);
}

size_t SimpleSolver::solve(const Graph in_graph, Graph& out_solutions)
{
    auto buf = in_graph.GetVerticesOrderedByDegree();
    m_p.clear();
    for (Vertex v : buf)
    {
        m_p.push_back(v);
    }
    min_color_num = m_p.size();
    SIZE = m_p.size();
    m_current_graph = in_graph;
    m_current_graph.GenerateNotNeighbours();
    RunGreedyHeuristic();
    std::vector<Vertex> init;
    init.clear();
    init.reserve(m_p.size());
    //findClique(init.data(), init.size(), buf.data(), buf.size());
    solveLP();
    checkClique(m_curr_clique, in_graph);
    out_solutions.SetVertexSet(m_curr_clique);
    return m_curr_clique.size();
}

void SimpleSolver::addConstraint(const std::vector<double>& constr_coef)
{
    size_t count = 0;
    for (Vertex v = 0; v < constr_coef.size(); ++v)
    {
        if (almostEqual(constr_coef[v], 1.))
        {
            count++;
            rmatval.push_back(1);
            rmatind.push_back(v);
        }
    }
    rmatbeg.push_back(m_prev_offset);
    m_prev_offset += count;
}

void SimpleSolver::solveLP() noexcept
{
    std::vector<Vertex> init;
    init.clear();
    init.reserve(m_p.size());
    generateConstraints();
    if (m_curr_clique.size() == min_color_num)
        return;
    m_obj.resize(m_p.size(), 1.);
    /* Declare and allocate space for the variables and arrays where we
    will store the optimization results including the status, objective
    value, variable values, dual values, row slacks and variable
    reduced costs. */

    int      solstat;
    double   objval;
    std::vector<double> solution;
    std::vector<double> pi;
    std::vector<double> slack;
    std::vector<double> dj;
    
    int           status = 0;
    int           i, j;
    int           cur_numrows, cur_numcols;
    /* Initialize the CPLEX environment */
    env = CPXopenCPLEX(&status);

    //For debug
    //CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON);
    //CPXsetintparam(env, CPXPARAM_Read_DataCheck, CPX_DATACHECK_WARN);
    //
    lp = CPXcreateprob(env, &status, "BnCProb");
    CPXchgobjsen(env, lp, CPX_MAX);

    CPXnewcols(env, lp, SIZE, m_obj.data(), m_lower_bounds.data(), m_upper_bounds.data(), NULL, NULL);

    senses.resize(NUMROWS);
    std::generate(senses.begin(), senses.end(), []() { return 'L'; });

    size_t NUMNZ = rmatval.size();
    CPXaddrows(env, lp, 0, NUMROWS, NUMNZ, m_main_c.data(), senses.data(), (const int*)rmatbeg.data(),
        (const int*)rmatind.data(), (const double*)rmatval.data(), NULL, NULL);

    status = CPXlpopt(env, lp);
    cur_numrows = CPXgetnumrows(env, lp);
    cur_numcols = CPXgetnumcols(env, lp);

    solution.resize(cur_numcols);
    slack.resize(cur_numrows);
    dj.resize(cur_numcols);
    pi.resize(cur_numrows);
    status = CPXgetx(env, lp, solution.data(), 0, cur_numcols - 1);
    status = CPXgetobjval(env, lp, &objval);
    if (objval <= m_curr_clique.size())
    {
        return;
    }

    Vertex for_branch;
    if (!getBranchVariable(solution, for_branch))
    {
        // all variables are integer, update clique
        m_curr_clique.clear();
        for (Vertex v : m_p)
        {
            if (almostEqual(solution[v], 1.))
            {
                m_curr_clique.push_back(v);
            }
        }
        return;
    }
    //branch = 1
    branchLP(solution, for_branch, 1.);

    //branch = 0
    branchLP(solution, for_branch, 0.);
}

void SimpleSolver::generateConstraints() noexcept
{
    m_lower_bounds.resize(m_p.size(), 0.);
    m_upper_bounds.resize(m_p.size(), 1.);
    generateColoringConstraints();
    generateIndependentSetsConstraints();
    generateEdgeConstraints();

    for each (auto& var in m_main_constraints)
    {
        NUMROWS++;
        addConstraint(var);
    }
    m_main_c.resize(m_main_constraints.size(), 1.);
}

void SimpleSolver::generateEdgeConstraints() noexcept
{
    for (Vertex i = 0; i < m_p.size(); ++i)
    {

        const auto& adjacent = m_current_graph.GetAdjacencyMatrixRow(i);
        for (Vertex j = i + 1; j < m_p.size(); ++j)
        {
            if (!adjacent[j])
            {
                if (already_added.find(std::make_pair(i, j)) != already_added.end() || already_added.find(std::make_pair(j, i)) != already_added.end())
                    continue;
                std::vector<double> constr_coef;
                constr_coef.resize(m_p.size(), 0.);
                constr_coef[i] = 1;
                constr_coef[j] = 1;
                m_main_constraints.push_back(std::move(constr_coef));
            }
        }
    }
}

void SimpleSolver::expandIndependentSet(Vertices& vert)
{
    std::set<Vertex> expand = m_current_graph.GetNotNeighbours(*(vert.begin()));
    for (Vertex v : vert)
    {
        expand = expand*m_current_graph.GetNotNeighbours(v);
    }
    while(!expand.empty())
    {
        Vertex v_max = *(expand.begin());
        size_t max_size = 0;
        for (const Vertex& v : expand)
        {
            const std::set<Vertex> & candidate = expand*m_current_graph.GetNotNeighbours(v);
            if (candidate.size() > max_size)
            {
                max_size = candidate.size();
                v_max = v;
            }
        }
        expand = expand*m_current_graph.GetNotNeighbours(v_max);
        vert.push_back(v_max);
    }    
}

void SimpleSolver::generateColoringConstraints() noexcept
{
    std::vector<GraphColoring> colorings;
    m_color_solver.solve(m_current_graph, colorings);
    for (auto& col : colorings)
    {
        //validateColoring(m_current_graph, col);
        size_t col_num = col.getColorCount();
        if (col_num < min_color_num)
            min_color_num = col_num;
        for (size_t color = 0; color < col_num; ++color)
        {
            auto & vert_by_col = col.getVerticesByColor(color);
            expandIndependentSet(vert_by_col);
            if (vert_by_col.size() <= 2)
                continue;
            std::vector<double> constr_coef;
            constr_coef.resize(m_p.size(), 0.);
            for (auto vert : vert_by_col)
            {
                constr_coef[vert] = 1;
                for (auto vert2 : vert_by_col)
                {
                    if (vert2 == vert)
                        continue;
                    if (already_added.find(std::make_pair(vert, vert2)) == already_added.end() && already_added.find(std::make_pair(vert2, vert)) == already_added.end())
                        already_added.insert(std::make_pair(vert, vert2));
                }
            }
            m_main_constraints.push_back(std::move(constr_coef));
        }
    }
}

void SimpleSolver::generateIndependentSetsConstraints() noexcept
{

    auto ordered_vert = m_current_graph.GetVerticesOrderedByDegree();
    std::reverse(ordered_vert.begin(), ordered_vert.end());
    std::set<Vertex> candidates(ordered_vert.begin(), ordered_vert.end());
    //std::copy(ordered_vert.begin(), ordered_vert.end(), candidates.begin());
    std::vector<std::set<Vertex>> greedy_ind_sets;
    while (!candidates.empty())
    {
        std::set<Vertex> expand = m_current_graph.GetNotNeighbours(*(ordered_vert.begin()))*candidates;
        std::set<Vertex> curr_set;
        curr_set.insert(*(ordered_vert.begin()));
        candidates.erase(candidates.find(*(ordered_vert.begin())));
        ordered_vert.erase(ordered_vert.begin());
        while (!expand.empty())
        {
            Vertex v_max = *expand.begin();
            std::set<Vertex> candidate;
            size_t max_size = 0;
            for (auto v : expand)
            {
                candidate = expand*m_current_graph.GetNotNeighbours(v);
                if (candidate.size() > max_size)
                {
                    max_size = candidate.size();
                    v_max = v;
                }
            }
            curr_set.insert(v_max);
            candidates.erase(candidates.find(v_max));
            ordered_vert.erase(std::find(ordered_vert.begin(), ordered_vert.end(), v_max));
            expand = expand*candidate*candidates;
        }
        if (curr_set.size() > 2)
            greedy_ind_sets.push_back(curr_set);
    }


    for (auto& greedy_set : greedy_ind_sets)
    {
        std::vector<double> constr_coef;
        constr_coef.resize(m_p.size(), 0.);
        for (auto vert : greedy_set)
        {
            constr_coef[vert] = 1;
            for (auto vert2 : greedy_set)
            {
                if (vert2 == vert)
                    continue;
                if (already_added.find(std::make_pair(vert, vert2)) == already_added.end() && already_added.find(std::make_pair(vert2, vert)) == already_added.end())
                    already_added.insert(std::make_pair(vert, vert2));
            }
        }
        m_main_constraints.push_back(std::move(constr_coef));
    }
}

bool SimpleSolver::getBranchVariable(const std::vector<double>& solution, Vertex & branch_var) noexcept
{
    double max_val = 0.0;
    bool found = false;
    branch_var = 0;
    for (Vertex& v : m_p)
    {
        if (!isIntegerVar(solution[v]) && solution[v] > max_val - std::numeric_limits<double>::epsilon())
        {
            branch_var = v;
            max_val = solution[v];
            found = true;
        }
    }
    return found;
}

void SimpleSolver::branchLP(std::vector<double>& solution, Vertex bvar, double bound) noexcept
{
    m_branches++;
    int           status = 0;
    double   objval;
    const char type_B[] = { 'B' };
    const char type_L[] = { 'L' };
    const char type_U[] = { 'U' };
    double normal_bounds[2] = { 0., 1. };
    status = CPXchgbds(env, lp, 1, (int*)&bvar, type_B, &bound);

    status = CPXlpopt(env, lp);
    CPXgetx(env, lp, solution.data(), 0, SIZE - 1);
    CPXgetobjval(env, lp, &objval);
    if (objval <= m_curr_clique.size())
    {
        return;
    }
    Vertex for_branch;
    if (!getBranchVariable(solution, for_branch))
    {
        // all variables are integer, update clique
        m_curr_clique.clear();
        for (Vertex v : m_p)
        {
            if (almostEqual(solution[v], 1.))
                m_curr_clique.push_back(v);
        }
        return;
    }
    //branch = 1
    branchLP(solution, for_branch, 1.);

    //branch = 0
    branchLP(solution, for_branch, 0.);


    status = CPXchgbds(env, lp, 1, (int*)&for_branch, type_U, normal_bounds+1);
    status = CPXchgbds(env, lp, 1, (int*)&for_branch, type_L, normal_bounds);

}

void SimpleSolver::RunGreedyHeuristic()
{
    std::vector<uint8_t> forbidden(m_current_graph.GetVertexCount(), 0);
    size_t max_times = m_current_graph.GetVertexCount() / 10;
    size_t sum = 0;
    for (size_t i = 0; i < max_times; ++i)
    {
        std::vector<Vertex> D;
        Vertices clique(0);
        for (Vertex v : m_p)
        {
            if (!forbidden[v])
                D.push_back(v);
        }
        if (m_curr_clique.size()> D.size() || forbidden.size()-sum < m_curr_clique.size())
            break;
        while (D.size())
        {
            Vertex k = D.front();
            D.erase(D.begin());
            clique.push_back(k);
            D = IntersectWithNeighbors(k, D.begin(), D.end());
        }
        if (clique.size() > m_curr_clique.size())
        {
            m_curr_clique = clique;
        }
        for (Vertex v : clique)
        {
            sum++;
            forbidden[v] = 1;
        }
    }
}

std::vector<Vertex> SimpleSolver::IntersectWithNeighbors(Vertex v, std::vector<Vertex>::const_iterator begin, std::vector<Vertex>::const_iterator end) noexcept
{
    std::vector<Vertex> out_P;
    out_P.reserve(std::distance(begin, end));
    const auto& adjacent = m_current_graph.GetAdjacencyMatrixRow(v);
    for (auto& i = begin; i<end;++i)
    {
        if (adjacent[*i])
        {
            out_P.push_back(*i);
        }
    }
    return out_P;
}


