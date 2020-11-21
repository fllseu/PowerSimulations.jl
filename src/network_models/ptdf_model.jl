function ptdf_networkflow(
    psi_container::PSIContainer,
    branches::IS.FlattenIteratorWrapper{B},
    buses::IS.FlattenIteratorWrapper{PSY.Bus},
    expression::Symbol,
    PTDF::PSY.PTDF,
) where {B <: PSY.Branch}
    time_steps = model_time_steps(psi_container)
    network_flow =
        add_cons_container!(psi_container, :network_flow, PTDF.axes[1], time_steps)

    constraint_val = JuMPConstraintArray(undef, time_steps)
    assign_constraint!(psi_container, "CopperPlateBalance", constraint_val)
    nodal_balance_expressions = psi_container.expressions[expression]

    remove_undef!(nodal_balance_expressions)

    _branches = sort!(
        collect(branches),
        by = x -> (PSY.get_number(PSY.get_arc(x).from), PSY.get_number(PSY.get_arc(x).to)),
    )

    typed_branches = IS.FlattenIteratorWrapper(
        PSY.ACBranch,
        Vector([_branches]),
    )
    flow_variables!(psi_container, StandardPTDFModel, typed_branches)


    for t in time_steps
        flow_variable = get_variable(psi_container, FLOW_ACTIVE_POWER, PSY.ACBranch)
        line_limit = (PTDF.data .* (nodal_balance_expressions.data[:, t]') ) * ones(length(PTDF.axes[2]),)

        network_flow[:, t] = JuMP.@constraint(
            psi_container.JuMPmodel,
            flow_variable[:, t].data .== line_limit
        )


        # The process is done in two separate loops to avoid modifying the nodal_balance_expressions
        # before making the flow constraints. If this two operations are done in the same loop
        # then the PTDF will multiply an expression that contains the flow variable.
        for br in branches
            name = PSY.get_name(br)
            from_number = PSY.get_number(PSY.get_arc(br).from)
            to_number = PSY.get_number(PSY.get_arc(br).to)
            flow_variable = get_variable(psi_container, FLOW_ACTIVE_POWER, PSY.ACBranch)
            add_to_expression!(
                nodal_balance_expressions,
                from_number,
                t,
                flow_variable[name, t],
                -1.0,
            )
            add_to_expression!(
                nodal_balance_expressions,
                to_number,
                t,
                flow_variable[name, t],
                1.0,
            )
        end

        constraint_val[t] = JuMP.@constraint(
            psi_container.JuMPmodel,
            sum(nodal_balance_expressions[:, t]) == 0
        )
    end
    return
end
