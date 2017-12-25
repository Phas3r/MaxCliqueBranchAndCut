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


double product(double* weights, const Vertices& vert)
{
    double sum = 0;
    for (const Vertex& v : vert)
    {
        sum += weights[v];
    }
    return sum;
}

double product(double* weights, const std::vector<int>& vert)
{
    double sum = 0;
    for (const Vertex& v : vert)
    {
        sum += weights[v];
    }
    return sum;
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
        still_unconstrained.insert(v);
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
    NUMROWS++;
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
    m_curr_non_zero.reserve(m_p.size());
    m_curr_most_violated.reserve(m_p.size());
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
    double   objval, old_objval;
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
    lp = CPXcreateprob(env, &status, "BnCProblem");
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


    char single_sence[] = { 'L' };
    double single_coef[] = { 1.0 };
    int single_matbeg[] = { 0, 0 };
    int single_matind[] = { 0 };
    status = CPXsetintparam(env, CPXPARAM_LPMethod, CPX_ALG_DUAL);
    while (findMostViolatedConstr(solution.data()))
    {
        status = CPXnewrows(env, lp, 1, NULL, NULL, NULL, NULL);
        status = CPXaddrows(env, lp, 0, 1, m_curr_most_violated.size(), m_main_c.data(), single_sence, single_matbeg,
            (const int*)m_curr_most_violated.data(), (const double*)m_upper_bounds.data(), NULL, NULL);


        NUMROWS++;
        cur_numrows = CPXgetnumrows(env, lp);

        status = CPXlpopt(env, lp);
        status = CPXsolution(env, lp, &solstat, &objval, solution.data(), 0, 0, 0);

        //status = CPXwriteprob(env, lp, "myprob.lp", NULL);
        if (objval <= m_curr_clique.size())
        {
            return;
        }
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
    generateIndependentSetsConstraints();
    generateColoringConstraints();
    //generateEdgeConstraints();

    for each (auto& var in m_main_constraints)
    {
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

void SimpleSolver::expandIndependentSet(Vertices& vert) noexcept
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

            for (auto v1 : vert_by_col)
            {
                auto g = m_current_graph.GetNotNeighbours(v1);
                for (auto v2 : vert_by_col)
                {
                    if (g.find(v2) == g.end() && v2 != v1)
                        std::cout << "KARAUL!";
                }
            }
            if (vert_by_col.size() <= 2)
                continue;

            std::vector<double> constr_coef;
            constr_coef.resize(m_p.size(), 0.);
            for (auto vert : vert_by_col)
            {
                still_unconstrained.erase(vert);
                constr_coef[vert] = 1;
                for (auto vert2 : vert_by_col)
                {
                    if (vert2 == vert)
                        continue;
                }
            }
            m_main_constraints.push_back(std::move(constr_coef));
        }
    }
    for (auto it = still_unconstrained.begin(); it != still_unconstrained.end(); ++it)
    {
        size_t max_size = 0;
        size_t max_ind = 0;
        for (size_t i = 0; i < colorings.size(); ++i)
        {
            if (colorings[i].getVerticesByColor(colorings[i].getVertexColor(*it)).size() > max_size)
            {
                max_size = colorings[i].getVerticesByColor(colorings[i].getVertexColor(*it)).size();
                max_ind = i;
            }
        }
        Vertices vert_by_col = colorings[max_ind].getVerticesByColor(colorings[max_ind].getVertexColor(*it));
        std::vector<double> constr_coef;
        constr_coef.resize(m_p.size(), 0.);
        for (auto vert : vert_by_col)
        {
            constr_coef[vert] = 1;
        }
        m_main_constraints.push_back(std::move(constr_coef));
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
            still_unconstrained.erase(vert);
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
    size_t added_in_branch = 0;
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


    while (findMostViolatedConstr(solution.data()))
    {

        status = CPXnewrows(env, lp, 1, NULL, NULL, NULL, NULL);
        status = CPXaddrows(env, lp, 0, 1, m_curr_most_violated.size(), m_main_c.data(), type_L, (const int*)rmatbeg.data(),
            (const int*)m_curr_most_violated.data(), (const double*)rmatval.data(), NULL, NULL);

        NUMROWS++;
        added_in_branch++;
        status = CPXlpopt(env, lp);
        status = CPXgetx(env, lp, solution.data(), 0, SIZE - 1);
        double pr = product(solution.data(), m_curr_most_violated);
        status = CPXgetobjval(env, lp, &objval);
        if (objval <= m_curr_clique.size())
        {
            status = CPXdelrows(env, lp, NUMROWS - added_in_branch, NUMROWS - 1);
            return;
        }
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
    status = CPXdelrows(env, lp, NUMROWS - added_in_branch, NUMROWS - 1);
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

bool SimpleSolver::findMostViolatedConstr(double* weights) noexcept
{
    m_curr_non_zero.clear();
    m_curr_most_violated.clear();
    for (size_t i = 0; i < SIZE; ++i)
    {
        if (!almostEqual(weights[i], 0.))
            m_curr_non_zero.push_back(i);
    }
    bool found = false;
    std::sort(m_curr_non_zero.begin(), m_curr_non_zero.end());
    GraphColoring rlf_coloring(m_curr_non_zero.size());
    size_t color_num = m_color_solver.solveRlf(m_current_graph, rlf_coloring, m_curr_non_zero);

    size_t needed_color = 0;
    double max_size = 1.0;
    for (size_t color = 0; color < color_num; ++color)
    {
        const Vertices& ind_set = rlf_coloring.getVerticesByColor(color);
        if (ind_set.size() < 2)
            continue;
        double pr = product(weights, ind_set);

        double sum = 0;
        for (const Vertex& v : ind_set)
        {
            sum += weights[m_curr_non_zero[v]];
        }
        if (sum > max_size)
        {
            max_size = pr;
            needed_color = color;
            found = true;
        }

        for (auto v1 : ind_set)
        {
            auto g = m_current_graph.GetNotNeighbours(m_curr_non_zero[v1]);
            for (auto v2 : ind_set)
            {
                if (g.find(m_curr_non_zero[v2]) == g.end() && m_curr_non_zero[v2] != m_curr_non_zero[v1])
                    std::cout << "KARAUL!";
            }
        }

    }
    if (found)
    {
        const Vertices& vert= rlf_coloring.getVerticesByColor(needed_color);
        for (auto v : vert)
        {
            m_curr_most_violated.push_back(m_curr_non_zero[v]);
        }
    }
    return found;
}
