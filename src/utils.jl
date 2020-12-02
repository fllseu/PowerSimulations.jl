const HASH_FILENAME = "check.sha256"

"""
Return a decoded JSON file.
"""
function read_json(filename::AbstractString)
    return open(filename) do io
        return JSON.parse(io)
    end
end

"""
Return the SHA 256 hash of a file.
"""
function compute_sha256(filename::AbstractString)
    return open(filename) do io
        return bytes2hex(SHA.sha256(io))
    end
end

"""
Return the key for the given value
"""
function find_key_with_value(d, value)
    for (k, v) in d
        v == value && return k
    end
    error("dict does not have value == $value")
end

function compute_file_hash(path::String, files::Vector{String})
    data = Dict("files" => [])
    for file in files
        file_path = joinpath(path, file)
        file_info = Dict("filename" => file_path, "hash" => compute_sha256(file_path))
        push!(data["files"], file_info)
    end

    open(joinpath(path, HASH_FILENAME), "w") do io
        write(io, JSON.json(data))
    end
end

function read_file_hashes(path)
    data = open(joinpath(path, HASH_FILENAME), "r") do io
        JSON.parse(io)
    end

    return data["files"]
end

function check_kwargs(input_kwargs, valid_set::Array{Symbol}, function_name::String)
    if isempty(input_kwargs)
        return
    else
        for (key, value) in input_kwargs
            if !(key in valid_set)
                throw(ArgumentError("keyword argument $(key) is not a valid input for $(function_name)"))
            end
        end
    end
    return
end

# writing a dictionary of dataframes to files

function write_data(vars_results::Dict, save_path::String; kwargs...)
    file_type = get(kwargs, :file_type, Arrow)
    if :duals in keys(kwargs)
        name = "dual_"
    elseif :params in keys(kwargs)
        name = "parameter_"
    else
        name = ""
    end
    if file_type == Arrow || file_type == CSV
        for (k, v) in vars_results
            file_path = joinpath(save_path, "$name$k.$(lowercase("$file_type"))")
            if isempty(vars_results[k])
                @debug "$name$k is empty, not writing $file_path"
            else
                file_type.write(file_path, vars_results[k])
            end
        end
    end
end

# writing a dictionary of dataframes to files and appending the time

function write_data(
    vars_results::Dict,
    time::DataFrames.DataFrame,
    save_path::AbstractString;
    kwargs...,
)
    file_type = get(kwargs, :file_type, Arrow)
    for (k, v) in vars_results
        var = DataFrames.DataFrame()
        if file_type == CSV && size(time, 1) == size(v, 1)
            var = hcat(time, v)
        else
            var = v
        end
        file_path = joinpath(save_path, "$(k).$(lowercase("$file_type"))")
        file_type.write(file_path, var)
    end
end

function write_data(
    data::DataFrames.DataFrame,
    save_path::AbstractString,
    file_name::String;
    kwargs...,
)
    if isfile(save_path)
        save_path = dirname(save_path)
    end
    file_type = get(kwargs, :file_type, Arrow)
    if file_type == Arrow || file_type == CSV
        file_path = joinpath(save_path, "$(file_name).$(lowercase("$file_type"))")
        file_type.write(file_path, data)
    end
    return
end

function write_optimizer_log(optimizer_log::Dict, save_path::AbstractString)
    JSON.write(joinpath(save_path, "optimizer_log.json"), JSON.json(optimizer_log))
end

function write_data(base_power::Float64, save_path::String)
    JSON.write(joinpath(save_path, "base_power.json"), JSON.json(base_power))
end

function _jump_value(input::JuMP.VariableRef)
    return JuMP.value(input)
end

function _jump_value(input::PJ.ParameterRef)
    return PJ.value(input)
end

function _jump_value(input::JuMP.ConstraintRef)
    return JuMP.dual(input)
end

function axis_array_to_dataframe(input_array::JuMP.Containers.DenseAxisArray{Float64})
    if length(axes(input_array)) == 1
        result = Vector{Float64}(undef, length(first(input_array.axes)))

        for t in input_array.axes[1]
            result[t] = input_array[t]
        end

        return DataFrames.DataFrame(var = result)

    elseif length(axes(input_array)) == 2
        result = Array{Float64, length(input_array.axes)}(
            undef,
            length(input_array.axes[2]),
            length(input_array.axes[1]),
        )
        names = Array{Symbol, 1}(undef, length(input_array.axes[1]))

        for t in input_array.axes[2], (ix, name) in enumerate(input_array.axes[1])
            result[t, ix] = input_array[name, t]
            names[ix] = Symbol(name)
        end

        return DataFrames.DataFrame(result, names)

    elseif length(axes(input_array)) == 3
        extra_dims = sum(length(axes(input_array)[2:(end - 1)]))
        extra_vars = [Symbol("S$(s)") for s in 1:extra_dims]
        result_df = DataFrames.DataFrame()
        names = vcat(extra_vars, Symbol.(axes(input_array)[1]))

        for i in input_array.axes[2]
            third_dim = collect(fill(i, size(input_array)[end]))
            result = Array{Float64, 2}(
                undef,
                length(last(input_array.axes)),
                length(first(input_array.axes)),
            )
            for t in last(input_array.axes),
                (ix, name) in enumerate(first(input_array.axes))

                result[t, ix] = input_array[name, i, t]
            end
            res = DataFrames.DataFrame(hcat(third_dim, result))
            result_df = vcat(result_df, res)
        end

        return DataFrames.rename!(result_df, names)

    else
        error("Dimension Number $(length(axes(input_array))) not Supported")
    end
end

function axis_array_to_dataframe(input_array::JuMP.Containers.DenseAxisArray{})
    if length(axes(input_array)) == 1
        result = Vector{Float64}(undef, length(first(input_array.axes)))
        for t in input_array.axes[1]
            result[t] = _jump_value(input_array[t])
        end

        return DataFrames.DataFrame(var = result)
    elseif length(axes(input_array)) == 2
        result = Array{Float64, length(input_array.axes)}(
            undef,
            length(input_array.axes[2]),
            length(input_array.axes[1]),
        )
        names = Array{Symbol, 1}(undef, length(input_array.axes[1]))

        for t in input_array.axes[2], (ix, name) in enumerate(input_array.axes[1])
            result[t, ix] = _jump_value(input_array[name, t])
            names[ix] = Symbol(name)
        end

        return DataFrames.DataFrame(result, names)

    elseif length(axes(input_array)) == 3
        extra_dims = sum(length(axes(input_array)[2:(end - 1)]))
        extra_vars = [Symbol("S$(s)") for s in 1:extra_dims]
        result_df = DataFrames.DataFrame()
        names = vcat(extra_vars, Symbol.(axes(input_array)[1]))

        for i in input_array.axes[2]
            third_dim = collect(fill(i, size(input_array)[end]))
            result = Array{Float64, 2}(
                undef,
                length(last(input_array.axes)),
                length(first(input_array.axes)),
            )
            for t in last(input_array.axes),
                (ix, name) in enumerate(first(input_array.axes))

                result[t, ix] = _jump_value(input_array[name, i, t])
            end
            res = DataFrames.DataFrame(hcat(third_dim, result))
            result_df = vcat(result_df, res)
        end
        return DataFrames.rename!(result_df, names)
    else
        @warn("Dimension Number $(length(axes(input_array))) not Supported, returning empty DataFrame")
        return DataFrames.DataFrame()
    end
end

function axis_array_to_dataframe(input_array::JuMP.Containers.SparseAxisArray)
    column_names = unique([(k[1], k[3]) for k in keys(input_array.data)])
    result_df = DataFrames.DataFrame()
    array_values = Vector{Vector{Float64}}(undef, length(column_names))
    for (ix, col) in enumerate(column_names)
        res = values(filter(v -> first(v)[[1, 3]] == col, input_array.data))
        array_values[ix] = PSI._jump_value.(res)
    end
    return DataFrames.DataFrame(array_values, Symbol.(column_names))
end

# this ensures that the time_stamp is not double shortened
function find_var_length(es::Dict, e_list::Array)
    return size(es[Symbol(splitext(e_list[1])[1])], 1)
end

function shorten_time_stamp(time::DataFrames.DataFrame)
    time = time[1:(size(time, 1) - 1), :]
    return time
end

""" Returns the correct container spec for the selected type of JuMP Model"""
function container_spec(m::M, axs...) where {M <: JuMP.AbstractModel}
    return JuMP.Containers.DenseAxisArray{JuMP.variable_type(m)}(undef, axs...)
end

""" Returns the correct container spec for the selected type of JuMP Model"""
function sparse_container_spec(m::M, axs...) where {M <: JuMP.AbstractModel}
    indexes = Base.Iterators.product(axs...)
    contents = Dict{eltype(indexes), Any}(indexes .=> 0)
    return JuMP.Containers.SparseAxisArray(contents)
end

function middle_rename(original::Symbol, split_char::String, addition::String)
    parts = split(String(original), split_char)
    parts[1] = parts[1] * "_" * addition
    return Symbol(join(parts, split_char))
end

"Replaces the string in `char` with the string`replacement`"
function replace_chars(s::String, char::String, replacement::String)
    return replace(s, Regex("[$char]") => replacement)
end

"Removes the string `char` from the original string"
function remove_chars(s::String, char::String)
    return replace_chars(s::String, char::String, "")
end

function get_available_components(::Type{T}, sys::PSY.System) where {T <: PSY.Component}
    return PSY.get_components(T, sys, x -> PSY.get_available(x))
end

function get_available_components(
    ::Type{PSY.RegulationDevice{T}},
    sys::PSY.System,
) where {T <: PSY.Component}
    return PSY.get_components(
        PSY.RegulationDevice{T},
        sys,
        x -> (PSY.get_available(x) && PSY.has_service(x, PSY.AGC)),
    )
end

"""
Load the complete arrow file into a DataFrame. Not optimized for memory use
"""
function read_arrow_file(file::AbstractString)
    return open(file, "r") do io
        DataFrames.DataFrame(Arrow.Table(io))
    end
end

function encode_symbol(::Type{T}, name1::AbstractString, name2::AbstractString) where {T}
    return Symbol(join((name1, name2, IS.strip_module_name(T)), PSI_NAME_DELIMITER))
end

function encode_symbol(
    ::Type{T},
    name1::AbstractString,
    name2::AbstractString,
) where {T <: PSY.Reserve}
    T_ = replace(IS.strip_module_name(T), "{" => "_")
    T_ = replace(T_, "}" => "")
    return Symbol(join((name1, name2, T_), PSI_NAME_DELIMITER))
end

function encode_symbol(::Type{T}, name1::Symbol, name2::Symbol) where {T}
    return encode_symbol(IS.strip_module_name(T), string(name1), string(name2))
end

function encode_symbol(::Type{T}, name::AbstractString) where {T}
    return Symbol(join((name, IS.strip_module_name(T)), PSI_NAME_DELIMITER))
end

function encode_symbol(::Type{T}, name::AbstractString) where {T <: PSY.Reserve}
    T_ = replace(IS.strip_module_name(T), "{" => "_")
    T_ = replace(T_, "}" => "")
    return Symbol(join((name, T_), PSI_NAME_DELIMITER))
end

function encode_symbol(::Type{T}, name::Symbol) where {T}
    return encode_symbol(T, string(name))
end

function encode_symbol(name::AbstractString)
    return Symbol(name)
end

function encode_symbol(name1::AbstractString, name2::AbstractString)
    return Symbol(join((name1, name2), PSI_NAME_DELIMITER))
end

function encode_symbol(name::Symbol)
    return name
end

function decode_symbol(name::Symbol)
    return split(String(name), PSI_NAME_DELIMITER)
end
