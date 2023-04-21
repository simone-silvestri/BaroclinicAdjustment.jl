using Oceananigans.Operators: ζ₃ᶠᶠᶜ
using Oceananigans.AbstractOperations: KernelFunctionOperation

using Oceananigans.Models: AbstractModel
using Oceananigans.Distributed

function standard_outputs!(simulation, output_prefix; overwrite_existing = true, 
                                                      checkpoint_time    = 100days,
                                                      snapshot_time      = 10days,
                                                      surface_time       = 1days)

    model = simulation.model
    grid  = model.grid

    u, v, w = model.velocities
    b = model.tracers.b

    output_fields = (; u, v, w, b)

    u2 = u^2
    v2 = v^2
    b2 = b^2
    w2 = w^2
    vb = v * b
    ub = u * b
    wb = w * b

    ζ  = KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid, u, v)
    ζ2 = ζ^2

    snapshot_fields = (; u, v, w, b, ζ, ζ2, u2, v2, w2, b2, ub, vb, wb)

    simulation.output_writers[:snapshots] = JLD2OutputWriter(model, snapshot_fields;
                                                                  schedule = TimeInterval(snapshot_time),
                                                                  filename = output_prefix * "_snapshots",
                                                                  overwrite_existing)

    simulation.output_writers[:surface_fields] = JLD2OutputWriter(model, output_fields;
                                                                  schedule = TimeInterval(surface_time),
                                                                  filename = output_prefix * "_surface",
                                                                  indices = (:, :, grid.Nz),
                                                                  overwrite_existing)

    simulation.output_writers[:checkpointer] = Checkpointer(model;
                                                            schedule = TimeInterval(checkpoint_time),
                                                            prefix = output_prefix * "_checkpoint",
                                                            overwrite_existing)

    return nothing
end

function checkpoint_outputs!(simulation, output_prefix; overwrite_existing = true, checkpoint_time = 100days)

    model = simulation.model

    simulation.output_writers[:checkpointer] = Checkpointer(model;
                                                            schedule = TimeInterval(checkpoint_time),
                                                            prefix = output_prefix * "_checkpoint",
                                                            overwrite_existing)

    return nothing
end

function reduced_outputs!(simulation, output_prefix; overwrite_existing = true, 
                                                     checkpoint_time    = 100days,
                                                     snapshot_time      = 1days,
                                                     surface_time       = 0.5days,
                                                     bottom_time        = 0.5days)

    model = simulation.model
    grid  = model.grid

    u, v, w = model.velocities
    b = model.tracers.b

    output_fields = (; u, v, w, b)

    simulation.output_writers[:snapshots] = JLD2OutputWriter(model, output_fields;
                                                                schedule = TimeInterval(snapshot_time),
                                                                filename = output_prefix * "_snapshots",
                                                                overwrite_existing)
   

    simulation.output_writers[:surface_fields] = JLD2OutputWriter(model, output_fields;
                                                                    schedule = TimeInterval(surface_time),
                                                                    filename = output_prefix * "_surface",
                                                                    indices = (:, :, grid.Nz),
                                                                    overwrite_existing)
                                                                                                                                
    simulation.output_writers[:bottom_fields] = JLD2OutputWriter(model, output_fields;
                                                                    schedule = TimeInterval(bottom_time),
                                                                    filename = output_prefix * "_bottom",
                                                                    indices = (:, :, 2),
                                                                    overwrite_existing)
    
    simulation.output_writers[:checkpointer] = Checkpointer(model;
                                                            schedule = TimeInterval(checkpoint_time),
                                                            prefix = output_prefix * "_checkpoint",
                                                            overwrite_existing)

end                                                 

