
import StaticArrays: SizedMatrix, SVector, MVector
import LinearAlgebra: cholesky, det, diag, tr, Diagonal


function set_sform!(header, sform::AbstractMatrix)
    T = typeof(header.srow_x)
    sform = SizedMatrix{4, 4}(sform)
    header.srow_x = T(sform[1, :])
    header.srow_y = T(sform[2, :])
    header.srow_z = T(sform[3, :])
    return header
end

function set_qform!(header::NiftiHeader, qform::AbstractMatrix)
    qform = SizedMatrix{4, 4}(qform)
    qoffset = qform[1:3, 4]
    quatern = mat2quat(qform[1:3, 1:3])
    pixdim = MVector(header.pixdim...)
    pixdim[1:4] = SVector(quatern[4:7])
    header.pixdim = typeof(header.pixdim)(pixdim)
    header.quatern = typeof(header.quatern)(quatern[1:3])
    header.qoffset = typeof(header.qoffset)(qoffset)
    return header
end

function set_qform!(header::NIfTI.NIfTI1Header, qform::AbstractMatrix)
    qform = SizedMatrix{4, 4}(qform)
    qoffset = qform[1:3, 4]
    quatern = mat2quat(qform[1:3, 1:3])
    pixdim = MVector(header.pixdim...)
    pixdim[1:4] = SVector(quatern[4:7])
    header.pixdim = typeof(header.pixdim)(pixdim)
    header.quatern_b = typeof(header.quatern_b)(quatern[1])
    header.quatern_c = typeof(header.quatern_c)(quatern[2])
    header.quatern_d = typeof(header.quatern_d)(quatern[3])
    header.qoffset_x = typeof(header.qoffset_x)(qoffset[1])
    header.qoffset_y = typeof(header.qoffset_y)(qoffset[2])
    header.qoffset_z = typeof(header.qoffset_z)(qoffset[3])
    return header
end


function mat2quat(M::AbstractMatrix)
    # extract rotation/scale
    flip = sign(det(M))
    if flip < 0
        M = copy(M)
        M[:, 3] *= -1
    end
    U = cholesky(M' * M).U
    M = M / U
    pixdim = diag(U)
    pixdim = (pixdim[1], pixdim[2], pixdim[3])
    # precompute values
    t = tr(M)
    d = diag(M)
    i = argmax(d)                 # index of largest diagonal entry
    m = maximum(abs.(d)) * 1E-9   # clip threshold
    if t > 0
        r = sqrt(abs(1 + t))
        s = 0.5 / clamp(r, m, Inf)
        w = 0.5 * r
        x = s * (M[3, 2] - M[2, 3])
        y = s * (M[1, 3] - M[3, 1])
        z = s * (M[2, 1] - M[1, 2])
    elseif i == 1
        r = sqrt(abs(1 + 2 * M[1, 1] - t))
        s = 0.5 / clamp(r, m, Inf)
        w = s * (M[3, 2] - M[2, 3])
        x = 0.5 * r
        y = s * (M[2, 1] + M[1, 2])
        z = s * (M[1, 3] + M[3, 1])
    elseif i == 2
        r = sqrt(abs(1 + 2 * M[2, 2] - t))
        s = 0.5 / clamp(r, m, Inf)
        w = s * (M[3, 1] - M[1, 3])
        x = s * (M[1, 2] + M[2, 1])
        y = 0.5 * r
        z = s * (M[3, 2] + M[2, 3])
    else  # i == 3
        r = sqrt(abs(1 + 2 * M[3, 3] - t))
        s = 0.5 / clamp(r, m, Inf)
        w = s * (M[2, 1] - M[1, 2])
        x = s * (M[1, 3] + M[3, 1])
        y = s * (M[3, 2] + M[2, 3])
        z = 0.5 * r
    end
    return (x, y, z, flip, pixdim...)
end
