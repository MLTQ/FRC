"""Latest saved fields and research status as a standalone scientific PNG."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from mpl_toolkits.mplot3d.art3d import Line3DCollection


def render(source, summary, data, target):
    """Plot hash-verified inputs supplied by current_state_display.build."""
    frame = summary['frames'][-1]
    with np.load(source/frame['arrays_file']) as raw:
        arrays = {k: raw[k] for k in raw.files}
    r, z = arrays['r'], arrays['z']
    current = np.linalg.norm(arrays['J'], axis=-1)/1000
    points, strengths = arrays['line_points'], arrays['line_strength_mT']
    scale = data['uniform_display_scale_per_m']
    points = points/scale*1000
    warm = LinearSegmentedColormap.from_list('current', [(0, '#210d24'), (.28, '#af264b'), (.62, '#f58442'), (1, '#fff5b9')])
    cool = LinearSegmentedColormap.from_list('magnetic', [(0, '#346db5'), (.45, '#39ccca'), (1, '#ffe2a0')])
    norms = [Normalize(0, data['jmax']/1000), Normalize(data['bmin_mT'], data['bmax_mT'])]
    fig = plt.figure(figsize=(16, 10.4), facecolor='#08101a')
    fig.text(.04, .947, 'VIRTUAL PLASMA BENCH  /  CURRENT RESEARCH CHECKPOINT', color='#8eb7bd', fontsize=10, family='monospace')
    fig.text(.039, .895, 'The latest saved plasma fields.', color='#e5f1ed', fontsize=30)
    fig.text(.04, .851, '28.33 µs · last accepted step 17     /     70 mm diameter × 240 mm     /     5 mT axial bias', color='#a4bcc9', fontsize=12)
    axes = [fig.add_axes([.025+i*.49, .29, .475, .48], projection='3d', computed_zorder=False) for i in range(2)]
    for i, ax in enumerate(axes):
        ax.set_facecolor('#08101a'); ax.set_axis_off()
        ax.set_xlim(-125, 125); ax.set_ylim(-38, 38); ax.set_zlim(-38, 38)
        ax.set_box_aspect((250, 76, 76), zoom=1.17); ax.view_init(elev=23, azim=-64)
        ax.set_proj_type('ortho')
        for end in (-120, 120):
            phi = np.linspace(0, 2*np.pi, 120)
            ax.plot(np.full_like(phi, end), 35*np.cos(phi), 35*np.sin(phi), color='#38505f', linewidth=.65, zorder=0)
        for phi in np.arange(6)*np.pi/3:
            ax.plot([-120, 120], [35*np.cos(phi)]*2, [35*np.sin(phi)]*2, color='#38505f', linewidth=.5, zorder=0)
        ax.text(-120, 0, -43, '−120 mm', color='#7896a8', fontsize=8)
        ax.text(96, 0, -43, '+120 mm  z', color='#7896a8', fontsize=8)
        fig.text(.055+i*.49, .802, ['01  CURRENT DENSITY', '02  MAGNETIC FIELD LINES'][i], color='#d4e7e5', fontsize=14)
        fig.text(.055+i*.49, .773, ['Cutaway surfaces · color = |J|', 'Full vector B · actual deflection · color = |B|'][i], color='#8da9b8', fontsize=10)
    ax = axes[0]
    ir = np.unique(np.r_[np.arange(0, len(r), 2), len(r)-1])
    iz = np.unique(np.r_[np.arange(0, len(z), 4), len(z)-1])
    rr, zz = np.meshgrid(r[ir]*1000, z[iz]*1000, indexing='ij')
    colors = warm(norms[0](current[np.ix_(ir, iz)]))
    # The +y,+vertical quarter is removed; both cut faces expose the scalar field.
    for phi in (0, np.pi/2):
        ax.plot_surface(zz, rr*np.cos(phi), rr*np.sin(phi), facecolors=colors,
                        rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False, zorder=3)
    phi = np.linspace(np.pi/2, 2*np.pi, 72)
    zz, pp = np.meshgrid(z[iz]*1000, phi, indexing='ij')
    outer = np.broadcast_to(current[-1, iz, None], zz.shape)
    ax.plot_surface(zz, r[-1]*1000*np.cos(pp), r[-1]*1000*np.sin(pp),
                    facecolors=warm(norms[0](outer)), rstride=1, cstride=1,
                    linewidth=0, antialiased=False, shade=False, alpha=.7, zorder=1)
    rr, pp = np.meshgrid(r[ir]*1000, phi, indexing='ij')
    for izend in (0, -1):
        colors = warm(norms[0](np.broadcast_to(current[ir, izend, None], rr.shape)))
        ax.plot_surface(np.full_like(rr, z[izend]*1000), rr*np.cos(pp), rr*np.sin(pp),
                        facecolors=colors, rstride=1, cstride=1, linewidth=0,
                        antialiased=False, shade=False, zorder=2)
    ax = axes[1]
    for k in range(points.shape[1]):
        p = points[:, k]
        line = Line3DCollection(np.stack((p[:-1], p[1:]), axis=1), cmap=cool, norm=norms[1], linewidths=1.1, alpha=.88)
        line.set_array((strengths[:-1, k]+strengths[1:, k])*.5); ax.add_collection3d(line)
    for i, cmap in enumerate((warm, cool)):
        ca = fig.add_axes([.057+i*.49, .285, .39, .009])
        bar = fig.colorbar(plt.cm.ScalarMappable(norm=norms[i], cmap=cmap), cax=ca, orientation='horizontal')
        bar.outline.set_visible(False); bar.ax.tick_params(labelsize=9, colors='#9cb7c6', length=0, pad=8)
        bar.set_label(['Current magnitude · kA/m² · shared linear scale', 'Magnetic magnitude · mT · shared linear scale'][i], color='#9cb7c6', size=10, labelpad=10)
    fig.text(.045, .201, f"Sampled peak |J|   {data['jmax']/1000:.1f} kA/m²", color='#d5e8e2', fontsize=14)
    fig.text(.535, .201, f"Sampled Bz minimum   +{frame['sampled_Bzmin_mT']:.4f} mT", color='#d5e8e2', fontsize=14)
    fig.text(.045, .145, 'Run stopped: step 18 failed to converge. Requested duration was 100 µs.', color='#e6c18b', fontsize=13)
    fig.text(.045, .108, 'Contact verification fixture · positive axial field · no FRC formation demonstrated.', color='#a5becb', fontsize=12)
    fig.text(.045, .073, 'Latest: contribution and local-current audits fully replayed; rapid growth remains unresolved.', color='#a5becb', fontsize=11)
    fig.text(.045, .035, 'AXISYMMETRIC DATA · PERIODIC ENDS · LINE SPACING / SURFACES ARE PRESENTATION CHOICES · 12 OCT 2026', color='#7896a8', fontsize=8, family='monospace')
    fig.savefig(target, dpi=180, facecolor=fig.get_facecolor()); plt.close(fig)
