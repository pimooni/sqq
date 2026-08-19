"""Embedded Tcl source used to generate SQQ VMD render scripts."""

from __future__ import annotations

SQQ_CAGE_TCL = r"""# SQQ annotated cage and guest renderer for VMD.
__SQQ_RENDER_MANIFEST__
if {[llength [info commands ::SQQ::dispose_renderer]] > 0} {
    catch {::SQQ::dispose_renderer 1}
} elseif {[namespace exists ::SQQ]} {
    # Compatibility cleanup for scripts generated before dispose_renderer.
    namespace eval ::SQQ {
        catch {trace remove variable ::vmd_frame write ::SQQ::frame_changed}
        catch {trace remove variable ::vmd_pick_atom write ::SQQ::pick_atom_changed}
        catch {trace remove variable ::vmd_pick_atom write ::SQQ::pick_atom_event}
        catch {trace remove variable ::vmd_pick_graphics write ::SQQ::pick_graphics_changed}
        foreach variable_name {frame_after_id pick_after_id renderer_after_id} {
            if {[info exists $variable_name] && [set $variable_name] ne ""} {
                catch {after cancel [set $variable_name]}
            }
        }
        catch {::SQQ::clear_graphics}
        catch {::SQQ::clear_representations}
    }
}
namespace eval ::SQQ {
    variable molid -1
    variable owns_molecule 0
    variable gro_path ""
    variable xtc_path ""
    variable membership_path ""
    variable active_families {cage}
    variable custom_show_active 0
    variable active_targets
    variable representation_names {}
    variable representation_keys {}
    variable pending_representation_keys {}
    variable representation_atoms
    variable representation_kind
    variable representation_color
    variable representation_material
    variable representation_radius_tier
    variable representation_name_by_key
    variable representation_registered
    variable representation_style_snapshot
    variable representation_style_template
    variable frame_after_id ""
    variable renderer_after_id ""
    variable renderer_state loading
    variable renderer_error ""
    variable ready_callbacks {}
    variable membership_frame_count 0
    variable membership_seen_frames
    variable parsed_frames 0
    variable displayed_graph_mode "__unset__"
    variable label_visible 0
    variable pick_mode off
    variable selected_cages {}
    variable selected_guest ""
    variable atom_label_count 0
    variable pick_after_id ""
    variable label_after_ids {}
    variable pick_cage_rep_name ""
    variable pick_guest_rep_name ""
    variable graphics_ids {}
    variable group_keys
    variable group_atoms
    variable graph_mode
    variable color_overrides
    variable known_objects
    variable object_aliases
    variable cage_types
    variable cage_ids
    variable cage_centers
    variable track_types
    variable graphics_targets
    variable guest_keys
    variable guest_atoms
    variable guest_types
    variable atom_guest
    variable component_atom_role
    variable component_resname_role
    variable components_initialized
    variable topology_indexes {}
    variable topology_resnames {}
    foreach name {group_keys group_atoms graph_mode color_overrides known_objects object_aliases cage_types cage_ids cage_centers track_types graphics_targets guest_keys guest_atoms guest_types atom_guest component_atom_role component_resname_role active_targets representation_atoms representation_kind representation_color representation_material representation_radius_tier representation_name_by_key representation_registered representation_style_snapshot representation_style_template membership_seen_frames components_initialized} {
        catch {array unset $name}
        array set $name {}
    }
    set active_targets(cage) [list [list cage *]]
}

proc ::SQQ::dispose_renderer {{delete_molecule 1}} {
    variable frame_after_id
    variable label_after_ids
    variable molid
    variable owns_molecule
    variable pick_after_id
    variable ready_callbacks
    variable renderer_after_id
    variable renderer_state
    set renderer_state disposing
    foreach variable_name {frame_after_id pick_after_id renderer_after_id} {
        if {[set $variable_name] ne ""} {
            catch {after cancel [set $variable_name]}
            set $variable_name ""
        }
    }
    foreach callback_id $label_after_ids { catch {after cancel $callback_id} }
    set label_after_ids {}
    catch {trace remove variable ::vmd_frame write ::SQQ::frame_changed}
    catch {trace remove variable ::vmd_pick_atom write ::SQQ::pick_atom_changed}
    catch {trace remove variable ::vmd_pick_atom write ::SQQ::pick_atom_event}
    catch {trace remove variable ::vmd_pick_graphics write ::SQQ::pick_graphics_changed}
    catch {::SQQ::clear_graphics}
    catch {::SQQ::clear_representations}
    if {$delete_molecule && $owns_molecule && $molid >= 0 &&
        ![catch {set molecule_ids [molinfo list]}] && $molid in $molecule_ids} {
        catch {mol delete $molid}
    }
    set molid -1
    set owns_molecule 0
    set ready_callbacks {}
}

proc ::SQQ::add_member {frame source key atom_index} {
    variable group_keys
    variable group_atoms
    variable known_objects
    set known_objects($source,$key) 1
    set group_key "$frame,$source"
    set atom_key "$frame,$source,$key"
    if {![info exists group_atoms($atom_key)]} {
        lappend group_keys($group_key) $key
        set group_atoms($atom_key) {}
    }
    lappend group_atoms($atom_key) $atom_index
}

proc ::SQQ::register_object {source key} {
    variable known_objects
    set known_objects($source,$key) 1
}

proc ::SQQ::register_cage_object {frame family cage_type object_id} {
    variable object_aliases
    variable cage_types
    set id_source "${family}-id"
    ::SQQ::register_object $family $cage_type
    ::SQQ::register_object $id_source $object_id
    set cage_types($frame,$object_id) $cage_type
    if {[regexp {_(t[1-9][0-9]*)$} $object_id -> track_id]} {
        variable track_types
        set track_source "${family}-track"
        ::SQQ::register_object $track_source $track_id
        set track_types($frame,$track_id) $cage_type
    }
    set compact [string map [list "^" "" "-" "" "_" ""] $cage_type]
    if {$compact ne $cage_type} {
        set object_aliases($family,$compact) $cage_type
        if {![string match "${cage_type}_*" $object_id]} { return }
        set suffix [string range $object_id [string length $cage_type] end]
        set object_aliases($id_source,${compact}${suffix}) $object_id
    }
}

proc ::SQQ::register_cage {frame cage_type object_id} {
    ::SQQ::register_cage_object $frame cage $cage_type $object_id
}

proc ::SQQ::register_guest_cage {frame cage_type object_id} {
    ::SQQ::register_cage_object $frame guest $cage_type $object_id
}

proc ::SQQ::set_cage_center {frame cage_id cage_type center} {
    variable cage_centers
    variable cage_ids
    set object_id [::SQQ::cage_object_id $cage_type $cage_id]
    ::SQQ::register_cage $frame $cage_type $object_id
    set key "$frame,$object_id"
    if {[info exists cage_centers($key)] && $cage_centers($key) ne $center} {
        error "Conflicting SQQ centers for cage $object_id in frame $frame"
    }
    set cage_centers($key) $center
    set cage_ids($key) $cage_id
}

proc ::SQQ::add_guest_group {frame identifier resname indexes} {
    variable atom_guest
    variable guest_atoms
    variable guest_keys
    variable guest_types
    set key "$frame,$identifier"
    if {[info exists guest_atoms($key)]} {
        error "Duplicate SQQ guest metadata for $identifier in frame $frame"
    }
    set indexes [lsort -integer -unique $indexes]
    if {[llength $indexes] == 0} {
        error "Empty SQQ guest metadata for $identifier in frame $frame"
    }
    set guest_atoms($key) $indexes
    set guest_types($key) $resname
    lappend guest_keys($frame) $identifier
    foreach atom_index $indexes {
        set atom_key "$frame,$atom_index"
        if {[info exists atom_guest($atom_key)] &&
            $atom_guest($atom_key) ne $identifier} {
            error "SQQ atom $atom_index belongs to multiple guests in frame $frame"
        }
        set atom_guest($atom_key) $identifier
    }
}

proc ::SQQ::normalize_component_role {value} {
    set role [string tolower [string trim $value]]
    if {$role ni {water guest additive environment other}} {
        error "Invalid SQQ component role '$value'"
    }
    return $role
}

proc ::SQQ::add_component_group {frame role resname indexes} {
    variable component_atom_role
    variable component_resname_role
    variable object_aliases
    set role [::SQQ::normalize_component_role $role]
    set resname [string trim $resname]
    if {$resname eq "" || $resname eq "-"} {
        error "SQQ component metadata requires a residue name"
    }
    ::SQQ::register_object component-role $role
    ::SQQ::register_object component-resname $resname
    set object_aliases(component-resname,[string tolower $resname]) $resname
    if {[info exists component_resname_role($resname)] &&
        $component_resname_role($resname) ne $role} {
        set component_resname_role($resname) other
    } else {
        set component_resname_role($resname) $role
    }
    foreach atom_index [lsort -integer -unique $indexes] {
        set atom_key "$frame,$atom_index"
        if {[info exists component_atom_role($atom_key)] &&
            $component_atom_role($atom_key) ne $role} {
            error "SQQ atom $atom_index has conflicting component roles in frame $frame"
        }
        set component_atom_role($atom_key) $role
        ::SQQ::add_member $frame component-role $role $atom_index
        ::SQQ::add_member $frame component-resname $resname $atom_index
    }
}

proc ::SQQ::deduplicate_frame_memberships {frame} {
    variable group_keys
    variable group_atoms
    foreach source {cage-id cage-track guest-id guest-track phase cluster domain component-role component-resname} {
        set group_key "$frame,$source"
        if {![info exists group_keys($group_key)]} { continue }
        foreach key $group_keys($group_key) {
            set atom_key "$frame,$source,$key"
            set group_atoms($atom_key) [lsort -integer -unique $group_atoms($atom_key)]
        }
    }
}

proc ::SQQ::numbered_id {prefix value} {
    if {[regexp {^[0-9]+$} $value]} {
        scan $value %d number
        return [format "%s_%05d" $prefix $number]
    }
    return [string tolower $value]
}

proc ::SQQ::cage_object_id {cage_type cage_id} {
    if {[regexp {^[0-9]+$} $cage_id]} {
        scan $cage_id %d number
        return [format "%s_%05d" $cage_type $number]
    }
    return "${cage_type}_$cage_id"
}

proc ::SQQ::cage_type_from_object_id {frame object_id} {
    variable cage_types
    variable track_types
    if {[info exists cage_types($frame,$object_id)]} {
        return $cage_types($frame,$object_id)
    }
    if {[info exists track_types($frame,$object_id)]} {
        return $track_types($frame,$object_id)
    }
    if {[regexp {^(.+)_([0-9]+)$} $object_id -> cage_type number]} {
        return $cage_type
    }
    return ""
}

proc ::SQQ::read_memberships {frame atom_index family payload} {
    if {$payload eq "" || $payload eq "-"} { return }
    foreach membership [split $payload ,] {
        set fields [split $membership :]
        if {[llength $fields] != 5} {
            error "Invalid SQQ $family membership in frame $frame: $membership"
        }
        lassign $fields cage_id cage_type phase domain_id cluster_id
        if {$cage_type eq "-"} { continue }
        set object_id [::SQQ::cage_object_id $cage_type $cage_id]
        if {$family eq "cage"} {
            ::SQQ::register_cage $frame $cage_type $object_id
            ::SQQ::add_member $frame cage-id $object_id $atom_index
            if {[regexp {^t[1-9][0-9]*$} $cage_id]} {
                ::SQQ::add_member $frame cage-track $cage_id $atom_index
            }
            if {$phase ne "-"} { ::SQQ::add_member $frame phase $phase $atom_index }
            if {$cluster_id ne "-"} {
                ::SQQ::add_member $frame cluster [::SQQ::numbered_id cluster $cluster_id] $atom_index
            }
            if {$domain_id ne "-"} {
                ::SQQ::add_member $frame domain [::SQQ::numbered_id domain $domain_id] $atom_index
            }
        } elseif {$family eq "guest"} {
            ::SQQ::register_guest_cage $frame $cage_type $object_id
            ::SQQ::add_member $frame guest-id $object_id $atom_index
            if {[regexp {^t[1-9][0-9]*$} $cage_id]} {
                ::SQQ::add_member $frame guest-track $cage_id $atom_index
            }
        } else {
            error "Invalid SQQ membership family '$family'"
        }
    }
}

proc ::SQQ::reset_membership_data {} {
    variable group_keys
    variable group_atoms
    variable graph_mode
    variable known_objects
    variable object_aliases
    variable cage_types
    variable cage_centers
    variable track_types
    variable guest_keys
    variable guest_atoms
    variable guest_types
    variable atom_guest
    variable component_atom_role
    variable component_resname_role
    variable components_initialized
    variable membership_frame_count
    variable membership_seen_frames
    variable topology_indexes
    variable topology_resnames
    foreach name {group_keys group_atoms graph_mode known_objects object_aliases cage_types cage_centers track_types guest_keys guest_atoms guest_types atom_guest component_atom_role component_resname_role components_initialized membership_seen_frames} {
        array unset $name
        array set $name {}
    }
    set membership_frame_count 0
    set topology_indexes {}
    set topology_resnames {}
}

proc ::SQQ::parse_membership_row {line line_number} {
    variable membership_frame_count
    variable membership_seen_frames
    set fields [split $line "\t"]
    if {[llength $fields] != 15} {
        error "Invalid SQQ membership TSV row $line_number"
    }
    lassign $fields record frame source_frame time_ps graph_value family cage_id cage_type phase domain_id cluster_id atom_indexes center_x center_y center_z
    if {![string is integer -strict $frame] || $frame < 0} {
        error "Invalid SQQ render frame at TSV row $line_number"
    }
    if {$record eq "F"} {
        variable graph_mode
        if {[info exists membership_seen_frames($frame)]} {
            error "Duplicate SQQ frame metadata for frame $frame"
        }
        set membership_seen_frames($frame) 1
        set graph_mode($frame) $graph_value
        if {$frame + 1 > $membership_frame_count} {
            set membership_frame_count [expr {$frame + 1}]
        }
    } elseif {$record eq "C"} {
        if {$family ne "cage" || $atom_indexes ne "-"} {
            error "Invalid SQQ cage-center record at TSV row $line_number"
        }
        foreach value [list $center_x $center_y $center_z] {
            if {![string is double -strict $value]} {
                error "Invalid SQQ cage center at TSV row $line_number"
            }
        }
        ::SQQ::set_cage_center $frame $cage_id $cage_type \
            [list $center_x $center_y $center_z]
    } elseif {$record eq "G"} {
        if {$family ne "guest" || $atom_indexes in {- ""}} {
            error "Invalid SQQ guest record at TSV row $line_number"
        }
        set indexes {}
        foreach atom_index [split $atom_indexes ,] {
            if {![string is integer -strict $atom_index] || $atom_index < 0} {
                error "Invalid SQQ guest atom index at TSV row $line_number"
            }
            lappend indexes $atom_index
        }
        ::SQQ::add_guest_group $frame $cage_id $cage_type $indexes
    } elseif {$record eq "P"} {
        if {$family ne "component" || $atom_indexes in {- ""}} {
            error "Invalid SQQ component record at TSV row $line_number"
        }
        # Render bundles have a fixed topology.  Frame-zero component metadata is
        # therefore canonical; later frames reuse it instead of duplicating a full
        # atom map in Tcl memory.
        if {$frame > 0} { return }
        set indexes {}
        foreach atom_index [split $atom_indexes ,] {
            if {![string is integer -strict $atom_index] || $atom_index < 0} {
                error "Invalid SQQ component atom index at TSV row $line_number"
            }
            lappend indexes $atom_index
        }
        ::SQQ::add_component_group $frame $cage_id $cage_type $indexes
    } elseif {$record eq "M"} {
        if {$family ni {cage guest}} {
            error "Invalid SQQ membership family at TSV row $line_number"
        }
        set payload "$cage_id:$cage_type:$phase:$domain_id:$cluster_id"
        if {$atom_indexes eq "-" || $atom_indexes eq ""} {
            error "Missing SQQ atom indexes at TSV row $line_number"
        }
        foreach atom_index [split $atom_indexes ,] {
            if {![string is integer -strict $atom_index] || $atom_index < 0} {
                error "Invalid SQQ atom index at TSV row $line_number"
            }
            ::SQQ::read_memberships $frame $atom_index $family $payload
        }
    } else {
        error "Invalid SQQ membership record at TSV row $line_number: $record"
    }
}

proc ::SQQ::read_membership_tsv {path} {
    ::SQQ::reset_membership_data
    variable membership_frame_count
    variable membership_seen_frames
    set handle [open $path r]
    fconfigure $handle -encoding ascii -translation auto
    if {[gets $handle header] < 0 || $header ne "record\trender_frame\tsource_frame\ttime_ps\tgraph_mode\tfamily\tcage_id\tcage_type\tphase\tdomain\tcluster\tatom_indices\tcenter_x_angstrom\tcenter_y_angstrom\tcenter_z_angstrom"} {
        close $handle
        error "Invalid SQQ membership TSV header: $path"
    }
    set line_number 1
    while {[gets $handle line] >= 0} {
        incr line_number
        if {$line eq ""} { continue }
        if {[catch {::SQQ::parse_membership_row $line $line_number} message]} {
            close $handle
            error $message
        }
    }
    close $handle
    for {set frame 0} {$frame < $membership_frame_count} {incr frame} {
        if {![info exists membership_seen_frames($frame)]} {
            error "Missing SQQ frame metadata for frame $frame"
        }
        ::SQQ::deduplicate_frame_memberships $frame
    }
    return $membership_frame_count
}

proc ::SQQ::renderer_fail {message} {
    variable renderer_error
    variable renderer_state
    set original_message $message
    ::SQQ::dispose_renderer 1
    set renderer_error $original_message
    set renderer_state failed
}

proc ::SQQ::require_renderer_ready {} {
    variable renderer_error
    variable renderer_state
    if {$renderer_state eq "ready"} { return }
    if {$renderer_state eq "failed"} {
        error "SQQ renderer failed: $renderer_error"
    }
    error "SQQ renderer is not ready"
}

proc ::SQQ::when_ready {script} {
    variable ready_callbacks
    variable renderer_error
    variable renderer_state
    if {$renderer_state eq "ready"} {
        uplevel #0 $script
    } elseif {$renderer_state eq "failed"} {
        error "SQQ renderer failed: $renderer_error"
    } else {
        lappend ready_callbacks $script
    }
}

proc ::SQQ::run_ready_callbacks {} {
    variable ready_callbacks
    set callbacks $ready_callbacks
    set ready_callbacks {}
    foreach script $callbacks {
        if {[catch {uplevel #0 $script} message options]} {
            puts stderr "SQQ deferred command failed: $message"
        }
    }
}

proc ::SQQ::initialize_component_frame {frame} {
    variable atom_guest
    variable component_atom_role
    variable components_initialized
    variable molid
    variable topology_indexes
    variable topology_resnames
    # Component identity is topology-level metadata.  Build it once for frame 0
    # and reuse the same atom indexes in every trajectory frame.
    if {$frame > 0} {
        if {![info exists components_initialized(0)]} {
            ::SQQ::initialize_component_frame 0
        }
        return
    }
    if {[info exists components_initialized(0)]} { return }
    if {[llength $topology_indexes] == 0} {
        if {[catch {
            set selection [atomselect $molid all frame 0]
            set topology_indexes [$selection get index]
            set topology_resnames [$selection get resname]
            $selection delete
        }]} {
            return
        }
    }
    if {[llength $topology_indexes] != [llength $topology_resnames]} {
        error "SQQ could not read component residue names from the render topology"
    }
    set has_metadata [expr {[llength [array names component_atom_role "${frame},*"]] > 0}]
    array set grouped {}
    foreach atom_index $topology_indexes resname $topology_resnames {
        set atom_key "$frame,$atom_index"
        if {[info exists component_atom_role($atom_key)]} {
            set role $component_atom_role($atom_key)
        } elseif {$frame > 0 && [info exists component_atom_role(0,$atom_index)]} {
            # A frame-zero component map is sufficient for a fixed topology.
            set role $component_atom_role(0,$atom_index)
        } elseif {$has_metadata} {
            set role other
        } elseif {[info exists atom_guest($atom_key)]} {
            # Compatibility fallback for compact bundles written before P records.
            set role guest
        } else {
            set role water
        }
        lappend grouped($role,$resname) $atom_index
    }
    foreach key [array names grouped] {
        lassign [split $key ,] role resname
        ::SQQ::add_component_group $frame $role $resname $grouped($key)
    }
    ::SQQ::deduplicate_frame_memberships $frame
    set components_initialized(0) 1
}

proc ::SQQ::load_renderer {} {
    variable gro_path
    variable membership_path
    variable molid
    variable owns_molecule
    variable parsed_frames
    variable xtc_path

    set parsed_frames [::SQQ::read_membership_tsv $membership_path]
    if {$parsed_frames == 0} {
        error "SQQ membership TSV contains no frames"
    }

    set molid [mol new $gro_path type gro waitfor all]
    set owns_molecule 1
    if {$molid < 0 || $molid ni [molinfo list]} {
        error "VMD discarded the SQQ topology while it was loading"
    }
    set topology_frames [molinfo $molid get numframes]
    if {$topology_frames != 1} {
        error "SQQ topology frame count mismatch: expected 1, VMD loaded $topology_frames"
    }
    if {$parsed_frames > 1} {
        mol addfile $xtc_path type xtc first 1 waitfor all molid $molid
    }
    if {$molid < 0 || $molid ni [molinfo list]} {
        error "VMD discarded the SQQ trajectory while it was loading"
    }
    set loaded [molinfo $molid get numframes]
    if {$loaded != $parsed_frames} {
        error "SQQ frame count mismatch: parsed $parsed_frames, VMD loaded $loaded"
    }
    ::SQQ::finish_renderer
}

proc ::SQQ::finish_renderer {} {
    variable displayed_graph_mode
    variable graph_mode
    variable molid
    variable parsed_frames
    variable ready_callbacks
    variable renderer_after_id
    variable renderer_state
    molinfo $molid set frame 0
    mol rename $molid __SQQ_MOLECULE_NAME__
    ::SQQ::initialize_component_frame 0
    ::SQQ::adopt_initial_representations
    display projection Orthographic
    # Keep startup output deterministic: graph-mode announcements begin only
    # when a later frame actually changes mode.
    set displayed_graph_mode [expr {[info exists graph_mode(0)] ?
        $graph_mode(0) : "unknown"}]
    ::SQQ::reset_show
    catch {color Display Background white}
    catch {trace remove variable ::vmd_frame write ::SQQ::frame_changed}
    trace add variable ::vmd_frame write ::SQQ::frame_changed
    catch {trace remove variable ::vmd_pick_atom write ::SQQ::pick_atom_event}
    trace add variable ::vmd_pick_atom write ::SQQ::pick_atom_event
    catch {trace remove variable ::vmd_pick_graphics write ::SQQ::pick_graphics_changed}
    trace add variable ::vmd_pick_graphics write ::SQQ::pick_graphics_changed
    set renderer_after_id ""
    set renderer_state ready
    puts "SQQ VMD Renderer: ready ($parsed_frames frame[expr {$parsed_frames == 1 ? "" : "s"}])"
    ::SQQ::startup_help
    if {[llength $ready_callbacks] > 0} {
        set renderer_after_id [after 1 [list ::SQQ::run_ready_callbacks]]
    }
}

proc ::SQQ::key_rank {family key} {
    if {$family in {cage guest}} {
        set order {512 51262 51263 51264 435663 51268}
        set position [lsearch -exact $order $key]
        if {$position >= 0} { return $position }
        return 100
    }
    if {$family eq "phase"} {
        set order {I II H B A U X}
        set position [lsearch -exact $order $key]
        if {$position >= 0} { return $position }
        return 100
    }
    return 0
}

proc ::SQQ::ordered_keys {family keys} {
    set decorated {}
    foreach key [lsort -unique $keys] {
        lappend decorated [list [::SQQ::key_rank $family $key] $key]
    }
    set output {}
    foreach item [lsort -integer -index 0 $decorated] { lappend output [lindex $item 1] }
    return $output
}

proc ::SQQ::standard_cage_rank {cage_type} {
    set order {512 51262 51263 51264 435663 51268}
    set position [lsearch -exact $order $cage_type]
    if {$position < 0} { return 0 }
    return [expr {$position + 1}]
}

proc ::SQQ::generic_cage_rank {cage_type} {
    set total_faces 0
    set ring_kinds 0
    set fields [regexp -all -inline {([0-9]+)\^([0-9]+)} $cage_type]
    for {set index 0} {$index < [llength $fields]} {incr index 3} {
        set count [lindex $fields [expr {$index + 2}]]
        if {[string is integer -strict $count]} {
            set total_faces [expr {$total_faces + $count}]
            incr ring_kinds
        }
    }
    return [list $total_faces $ring_kinds]
}

proc ::SQQ::object_render_key {frame object_id color_priority color_id explicit} {
    set cage_type [::SQQ::cage_type_from_object_id $frame $object_id]
    set exact [expr {$explicit || $color_priority == 3}]
    set standard_rank [::SQQ::standard_cage_rank $cage_type]
    if {$standard_rank > 0} {
        set standard 1
        set primary $standard_rank
        set secondary 0
    } else {
        set standard 0
        lassign [::SQQ::generic_cage_rank $cage_type] primary secondary
    }
    set exact_id [expr {$exact ? $object_id : ""}]
    return [list $exact $standard $primary $secondary $cage_type $exact_id $color_id]
}

proc ::SQQ::compare_object_render_keys {left right} {
    foreach index {0 1 2 3} {
        set left_value [lindex $left $index]
        set right_value [lindex $right $index]
        if {$left_value != $right_value} {
            return [expr {$left_value < $right_value ? -1 : 1}]
        }
    }
    foreach index {4 5} {
        set comparison [string compare [lindex $left $index] [lindex $right $index]]
        if {$comparison != 0} { return $comparison }
    }
    set left_color [lindex $left 6]
    set right_color [lindex $right 6]
    return [expr {$left_color < $right_color ? -1 : ($left_color > $right_color)}]
}

proc ::SQQ::cage_radius_tier {render_key} {
    if {[lindex $render_key 0]} { return 7 }
    if {[lindex $render_key 1]} { return [lindex $render_key 2] }
    return 0
}

proc ::SQQ::cage_layer_radius {tier tiers} {
    set count [llength $tiers]
    if {$count <= 1} { return 0.125 }
    set index [lsearch -exact $tiers $tier]
    return [format "%.3f" [expr {0.125 + 0.005 * $index / double($count - 1)}]]
}

proc ::SQQ::active_cage_radius_tiers {} {
    variable active_families
    variable active_targets
    variable color_overrides
    variable known_objects
    if {"cage" ni $active_families || ![info exists active_targets(cage)]} {
        return {}
    }
    array set tiers {}
    foreach target $active_targets(cage) {
        lassign $target source key
        if {$source in {cage-id cage-track}} {
            set tiers(7) 1
        } elseif {$source eq "cage" && $key eq "*"} {
            foreach name [array names known_objects "cage,*"] {
                set cage_type [string range $name 5 end]
                set rank [::SQQ::standard_cage_rank $cage_type]
                set tiers([expr {$rank > 0 ? $rank : 0}]) 1
            }
        } elseif {$source eq "cage"} {
            set rank [::SQQ::standard_cage_rank $key]
            set tiers([expr {$rank > 0 ? $rank : 0}]) 1
        }
    }
    foreach name [array names color_overrides] {
        set source [lindex [split $name ,] 0]
        if {$source in {cage-id cage-track}} { set tiers(7) 1 }
    }
    return [lsort -integer [array names tiers]]
}

proc ::SQQ::stable_color {key} {
    set palette {0 1 7 3 11 10 4 9 5 6}
    set hash 0
    foreach character [split $key ""] {
        scan $character %c code
        set hash [expr {(($hash * 33) + $code) & 0x7fffffff}]
    }
    return [lindex $palette [expr {$hash % [llength $palette]}]]
}

proc ::SQQ::color_id {family key} {
    if {$family in {cage guest}} {
        switch -- $key {
            512 { return 7 }
            51262 { return 0 }
            51263 { return 1 }
            51264 { return 3 }
            51268 { return 11 }
            435663 { return 10 }
            default { return 2 }
        }
    }
    if {$family eq "phase"} {
        switch -- $key {
            I { return 1 }
            II { return 0 }
            H { return 7 }
            B { return 3 }
            A { return 11 }
            U { return 2 }
            X { return 4 }
            default { return 2 }
        }
    }
    if {$family eq "component"} {
        switch -- [string tolower $key] {
            water { return 0 }
            guest { return 3 }
            additive { return 4 }
            environment { return 2 }
            other { return 6 }
        }
    }
    return [::SQQ::stable_color $key]
}

proc ::SQQ::source_family {source} {
    switch -- $source {
        cage - cage-id - cage-track { return cage }
        guest - guest-id - guest-track { return guest }
        phase - cluster - domain { return $source }
        component - component-role - component-resname { return component }
    }
    error "Unknown SQQ source '$source'"
}

proc ::SQQ::default_color_id {frame source key} {
    variable component_resname_role
    set family [::SQQ::source_family $source]
    if {$source in {cage-id cage-track guest-id guest-track}} {
        set cage_type [::SQQ::cage_type_from_object_id $frame $key]
        if {$cage_type ne ""} { return [::SQQ::color_id $family $cage_type] }
        return 2
    }
    if {$source eq "component-resname" && [info exists component_resname_role($key)]} {
        return [::SQQ::color_id component $component_resname_role($key)]
    }
    return [::SQQ::color_id $family $key]
}

proc ::SQQ::effective_color {frame source key} {
    variable color_overrides
    set exact "$source,$key"
    if {[info exists color_overrides($exact)]} {
        if {$color_overrides($exact) eq "default"} {
            return [list [::SQQ::default_color_id $frame $source $key] 3]
        }
        return [list $color_overrides($exact) 3]
    }
    set family [::SQQ::source_family $source]
    if {$source in {cage-id cage-track guest-id guest-track}} {
        set cage_type [::SQQ::cage_type_from_object_id $frame $key]
        if {$cage_type ne ""} {
            set type_key "$family,$cage_type"
            if {[info exists color_overrides($type_key)]} {
                if {$color_overrides($type_key) eq "default"} {
                    return [list [::SQQ::default_color_id $frame $source $key] 2]
                }
                return [list $color_overrides($type_key) 2]
            }
        }
    }
    set category_key "$family,*"
    if {[info exists color_overrides($category_key)]} {
        return [list $color_overrides($category_key) 1]
    }
    return [list [::SQQ::default_color_id $frame $source $key] 0]
}

proc ::SQQ::phase_key {value} {
    set aliases [dict create si I i I sii II ii II sh H h H boundary B b B ambiguous A a A unclassified U u U isolated X x X]
    set key [string tolower $value]
    if {[dict exists $aliases $key]} { return [dict get $aliases $key] }
    return ""
}

proc ::SQQ::normalize_family {value} {
    set family [string tolower [string trim $value]]
    if {$family ni {cage guest phase cluster domain component}} {
        error "SQQ family must be cage, guest, phase, cluster, domain, or component"
    }
    return $family
}

proc ::SQQ::is_family_token {value} {
    return [expr {[string tolower [string trim $value]] in {cage guest phase cluster domain component}}]
}

proc ::SQQ::parse_target {family value} {
    variable known_objects
    variable object_aliases
    set family [::SQQ::normalize_family $family]
    set token [string trim $value]
    if {[string equal -nocase $token all]} {
        if {$family eq "component"} { return [list component-role *] }
        return [list $family *]
    }
    if {$family in {cage guest}} {
        set id_source "${family}-id"
        set track_source "${family}-track"
        if {[regexp {^t[1-9][0-9]*$} $token] &&
            [info exists known_objects($track_source,$token)]} {
            return [list $track_source $token]
        }
        foreach source [list $id_source $family] {
            if {[info exists known_objects($source,$token)]} { return [list $source $token] }
            if {[info exists object_aliases($source,$token)]} {
                return [list $source $object_aliases($source,$token)]
            }
        }
        if {[regexp {^(.+)_([0-9]+)$} $token -> cage_type number]} {
            return [list $id_source [::SQQ::cage_object_id $cage_type $number]]
        }
        if {[regexp {^[0-9]+$} $token]} { return [list $family $token] }
    } elseif {$family eq "phase"} {
        set key [::SQQ::phase_key $token]
        if {$key ne ""} { return [list phase $key] }
    } elseif {$family eq "cluster"} {
        if {[regexp -nocase {^cluster_([0-9]+)$} $token -> number]} {
            return [list cluster [::SQQ::numbered_id cluster $number]]
        }
        if {[regexp {^[0-9]+$} $token]} {
            return [list cluster [::SQQ::numbered_id cluster $token]]
        }
    } elseif {$family eq "domain"} {
        if {[regexp -nocase {^domain_([0-9]+)$} $token -> number]} {
            return [list domain [::SQQ::numbered_id domain $number]]
        }
        if {[regexp {^[0-9]+$} $token]} {
            return [list domain [::SQQ::numbered_id domain $token]]
        }
    } elseif {$family eq "component"} {
        set role [string tolower $token]
        if {$role in {water guest additive environment other} &&
            [info exists known_objects(component-role,$role)]} {
            return [list component-role $role]
        }
        if {[info exists known_objects(component-resname,$token)]} {
            return [list component-resname $token]
        }
        set alias [string tolower $token]
        if {[info exists object_aliases(component-resname,$alias)]} {
            return [list component-resname $object_aliases(component-resname,$alias)]
        }
    }
    error "Unknown SQQ $family target '$value'"
}

proc ::SQQ::require_known_target {family target} {
    variable known_objects
    lassign $target source key
    if {$key eq "*"} { return }
    if {![info exists known_objects($source,$key)]} {
        error "SQQ $family target '[::SQQ::target_label $family $target]' does not exist in the loaded trajectory"
    }
}

proc ::SQQ::phase_label {key} {
    switch -- $key {
        I { return sI }
        II { return sII }
        H { return sH }
        B { return boundary }
        A { return ambiguous }
        U { return unclassified }
        X { return isolated }
    }
    return $key
}

proc ::SQQ::target_label {family target} {
    lassign $target source key
    if {$key eq "*"} { return all }
    if {$family eq "phase"} { return [::SQQ::phase_label $key] }
    return $key
}

proc ::SQQ::parse_targets {family values} {
    if {[llength $values] == 0} { error "At least one SQQ target is required" }
    set targets {}
    array set seen {}
    foreach value $values {
        set target [::SQQ::parse_target $family $value]
        ::SQQ::require_known_target $family $target
        if {![info exists seen($target)]} {
            set seen($target) 1
            lappend targets $target
        }
    }
    if {[llength $targets] > 1} {
        foreach target $targets {
            if {[lindex $target 1] eq "*"} { error "The all target must be used alone" }
        }
    }
    return $targets
}

proc ::SQQ::parse_show_groups {values} {
    if {[llength $values] < 2} {
        error "Usage: sqq show <family> <target...> ?<family> <target...> ...?"
    }
    set groups {}
    set family ""
    set raw_targets {}
    foreach value $values {
        if {[::SQQ::is_family_token $value]} {
            if {$family ne ""} {
                if {[llength $raw_targets] == 0} {
                    error "SQQ show family '$family' requires at least one target"
                }
                lappend groups [list $family [::SQQ::parse_targets $family $raw_targets]]
            }
            set family [::SQQ::normalize_family $value]
            set raw_targets {}
        } else {
            if {$family eq ""} {
                error "SQQ show must begin with cage, guest, phase, cluster, domain, or component"
            }
            lappend raw_targets $value
        }
    }
    if {$family eq "" || [llength $raw_targets] == 0} {
        if {$family eq ""} {
            error "SQQ show must begin with cage, guest, phase, cluster, domain, or component"
        }
        error "SQQ show family '$family' requires at least one target"
    }
    lappend groups [list $family [::SQQ::parse_targets $family $raw_targets]]
    return $groups
}

proc ::SQQ::merge_show_group {family targets} {
    variable active_families
    variable active_targets
    if {$family ni $active_families} {
        lappend active_families $family
        set active_targets($family) {}
    }
    set includes_all 0
    foreach target $targets {
        if {[lindex $target 1] eq "*"} {
            set includes_all 1
            break
        }
    }
    set already_all 0
    foreach target $active_targets($family) {
        if {[lindex $target 1] eq "*"} {
            set already_all 1
            break
        }
    }
    if {$includes_all} {
        set active_targets($family) $targets
    } elseif {!$already_all} {
        foreach target $targets {
            if {$target ni $active_targets($family)} {
                lappend active_targets($family) $target
            }
        }
    }
}

proc ::SQQ::set_show {values args} {
    variable active_families
    variable active_targets
    variable custom_show_active
    if {[llength $args] > 0} {
        set combined [list $values]
        foreach argument $args {
            foreach value $argument { lappend combined $value }
        }
        set values $combined
    }
    set groups [::SQQ::parse_show_groups $values]
    if {!$custom_show_active} {
        set active_families {}
        array unset active_targets
        array set active_targets {}
        set custom_show_active 1
    }
    foreach group $groups {
        lassign $group family targets
        ::SQQ::merge_show_group $family $targets
    }
    ::SQQ::render_current 1
}

proc ::SQQ::reset_show {{announce 0}} {
    ::SQQ::cancel_pending_pick
    variable active_families
    variable active_targets
    variable color_overrides
    variable custom_show_active
    variable label_visible
    variable pick_mode
    variable selected_cages
    variable selected_guest
    set active_families {cage}
    array unset active_targets
    array set active_targets {}
    set active_targets(cage) [list [list cage *]]
    array unset color_overrides
    array set color_overrides {}
    set custom_show_active 0
    set label_visible 0
    set pick_mode off
    set selected_cages {}
    set selected_guest ""
    ::SQQ::render_current
    if {$announce} { puts "SQQ clear: restored source-time cage view" }
}

proc ::SQQ::set_label {values} {
    variable label_visible
    set count [llength $values]
    if {$count == 0} {
        set label_visible [expr {!$label_visible}]
    } elseif {$count == 1} {
        set value [string tolower [string trim [lindex $values 0]]]
        if {$value ni {on off}} { error "Usage: sqq show label ?on|off?" }
        set label_visible [expr {$value eq "on"}]
    } else {
        error "Usage: sqq show label ?on|off?"
    }
    ::SQQ::update_current
    puts "SQQ label: [expr {$label_visible ? "on" : "off"}]"
}

proc ::SQQ::set_pick_mode {value} {
    ::SQQ::cancel_pending_pick
    variable atom_label_count
    variable pick_mode
    variable selected_cages
    variable selected_guest
    set mode [string tolower [string trim $value]]
    if {$mode ni {center guest off}} {
        error "Usage: sqq pick center|guest|off"
    }
    set pick_mode $mode
    set selected_cages {}
    set selected_guest ""
    ::SQQ::render_current
    if {$mode eq "off"} {
        puts "SQQ pick: off; restored opaque objects"
    } else {
        if {![catch {label list Atoms} atom_labels]} {
            set atom_label_count [llength $atom_labels]
        } else {
            set atom_label_count 0
        }
        set native_mode [expr {$mode eq "center" ? "labelatom" : "pick"}]
        if {[catch {mouse mode $native_mode} message]} {
            set pick_mode off
            ::SQQ::render_current
            error "SQQ could not enable VMD picking: $message"
        }
        puts "SQQ pick: $mode; VMD $native_mode mode enabled; objects are transparent until selected"
    }
}

proc ::SQQ::base_material {} {
    variable pick_mode
    return [expr {$pick_mode eq "off" ? "Opaque" : "Transparent"}]
}

proc ::SQQ::color_value {value} {
    if {[string equal -nocase $value default]} { return default }
    set names [colorinfo colors]
    if {[string is integer -strict $value]} {
        scan $value %d color_id
        if {$color_id < 0 || $color_id >= [llength $names]} {
            error "VMD ColorID must be between 0 and [expr {[llength $names] - 1}]"
        }
        return $color_id
    }
    set color_id [lsearch -nocase -exact $names $value]
    if {$color_id < 0} { error "Unknown VMD color '$value'" }
    return $color_id
}

proc ::SQQ::clear_family_colors {family} {
    variable color_overrides
    foreach name [array names color_overrides] {
        set source [lindex [split $name ,] 0]
        if {[::SQQ::source_family $source] eq $family} { unset color_overrides($name) }
    }
}

proc ::SQQ::set_colors {family values color} {
    variable color_overrides
    set family [::SQQ::normalize_family $family]
    set targets [::SQQ::parse_targets $family $values]
    set value [::SQQ::color_value $color]
    foreach target $targets {
        lassign $target source key
        if {$key eq "*"} {
            ::SQQ::clear_family_colors $family
            if {$value ne "default"} { set color_overrides($family,*) $value }
        } else {
            set color_overrides($source,$key) $value
        }
    }
    ::SQQ::render_current
    set labels {}
    foreach target $targets { lappend labels [::SQQ::target_label $family $target] }
    puts "SQQ color $family: [join $labels { }] -> $color"
}

proc ::SQQ::track_representation {rep_index} {
    variable molid
    variable representation_names
    if {[catch {mol repname $molid $rep_index} name]} { return }
    if {$name ni $representation_names} { lappend representation_names $name }
}

proc ::SQQ::adopt_initial_representations {} {
    variable molid
    set count [molinfo $molid get numreps]
    for {set rep 0} {$rep < $count} {incr rep} { ::SQQ::track_representation $rep }
}

proc ::SQQ::clear_representations {} {
    variable molid
    variable pick_cage_rep_name
    variable pick_guest_rep_name
    variable representation_atoms
    variable representation_color
    variable representation_keys
    variable representation_kind
    variable representation_material
    variable representation_name_by_key
    variable representation_names
    variable representation_registered
    variable representation_radius_tier
    variable representation_style_snapshot
    variable representation_style_template
    variable pending_representation_keys
    set indexes {}
    if {$molid >= 0 && ![catch {set molecule_ids [molinfo list]}] &&
        $molid in $molecule_ids} {
        foreach name $representation_names {
            if {[catch {mol repindex $molid $name} rep]} { continue }
            if {[string is integer -strict $rep] && $rep >= 0} { lappend indexes $rep }
        }
        foreach rep [lsort -integer -decreasing -unique $indexes] {
            catch {mol delrep $rep $molid}
        }
    }
    set representation_names {}
    set representation_keys {}
    set pending_representation_keys {}
    foreach name {representation_atoms representation_color representation_kind representation_material representation_name_by_key representation_radius_tier representation_registered representation_style_snapshot representation_style_template} {
        array unset $name
        array set $name {}
    }
    set pick_cage_rep_name ""
    set pick_guest_rep_name ""
}

proc ::SQQ::expanded_targets {frame family targets} {
    variable group_keys
    set expanded {}
    array set seen {}
    foreach target $targets {
        lassign $target source key
        if {$family in {cage guest} && $source eq $family} {
            set id_source "${family}-id"
            set group_key "$frame,$id_source"
            if {![info exists group_keys($group_key)]} { continue }
            foreach object_id [lsort -unique $group_keys($group_key)] {
                if {$key eq "*" || [::SQQ::cage_type_from_object_id $frame $object_id] eq $key} {
                    set item [list $id_source $object_id]
                    if {![info exists seen($item)]} {
                        set seen($item) 1
                        lappend expanded $item
                    }
                }
            }
        } elseif {$key eq "*"} {
            set group_key "$frame,$source"
            if {$family eq "component" && ![info exists group_keys($group_key)]} {
                set group_key "0,$source"
            }
            if {![info exists group_keys($group_key)]} { continue }
            foreach object_key [::SQQ::ordered_keys $family $group_keys($group_key)] {
                set item [list $source $object_key]
                if {![info exists seen($item)]} {
                    set seen($item) 1
                    lappend expanded $item
                }
            }
        } else {
            set item [list $source $key]
            if {![info exists seen($item)]} {
                set seen($item) 1
                lappend expanded $item
            }
        }
    }
    return $expanded
}

proc ::SQQ::group_atom_key {frame source key} {
    variable group_atoms
    set atom_key "$frame,$source,$key"
    if {[info exists group_atoms($atom_key)]} { return $atom_key }
    if {[::SQQ::source_family $source] eq "component"} {
        set topology_key "0,$source,$key"
        if {[info exists group_atoms($topology_key)]} { return $topology_key }
    }
    return $atom_key
}

proc ::SQQ::compare_render_keys {left right} {
    lassign $left left_priority left_color
    lassign $right right_priority right_color
    if {$left_priority != $right_priority} {
        return [expr {$left_priority < $right_priority ? -1 : 1}]
    }
    return [expr {$left_color < $right_color ? -1 : ($left_color > $right_color)}]
}

proc ::SQQ::representation_style_class {rep_key kind radius_tier} {
    set scope [lindex $rep_key 0]
    set family [lindex $rep_key 1]
    set semantic_type $kind
    if {$family in {cage guest}} {
        set layer_key [lindex $rep_key 2]
        if {[llength $layer_key] > 4} { set semantic_type [lindex $layer_key 4] }
    }
    return [list $scope $family $semantic_type $radius_tier]
}

proc ::SQQ::current_representation_style {rep} {
    variable molid
    if {[catch {
        set values [molinfo $molid get [list \
            [list rep $rep] [list color $rep] [list material $rep]]]
    }]} {
        return {}
    }
    if {[llength $values] != 3} { return {} }
    return $values
}

proc ::SQQ::refresh_representation_style_templates {} {
    variable molid
    variable representation_keys
    variable representation_kind
    variable representation_name_by_key
    variable representation_radius_tier
    variable representation_style_snapshot
    variable representation_style_template
    foreach rep_key $representation_keys {
        if {![info exists representation_name_by_key($rep_key)]} { continue }
        set rep_name $representation_name_by_key($rep_key)
        if {[catch {mol repindex $molid $rep_name} rep] || $rep < 0} { continue }
        set style [::SQQ::current_representation_style $rep]
        if {[llength $style] != 3} { continue }
        if {![info exists representation_style_snapshot($rep_key)]} {
            set representation_style_snapshot($rep_key) $style
            continue
        }
        if {$style ne $representation_style_snapshot($rep_key)} {
            set style_class [::SQQ::representation_style_class $rep_key \
                $representation_kind($rep_key) $representation_radius_tier($rep_key)]
            set representation_style_template($style_class) $style
            set representation_style_snapshot($rep_key) $style
        }
    }
}

proc ::SQQ::register_stable_representation {frame rep_key kind indexes color_id material {radius_tier ""}} {
    variable representation_atoms
    variable representation_color
    variable representation_keys
    variable representation_kind
    variable representation_material
    variable representation_registered
    variable representation_radius_tier
    variable pending_representation_keys
    set indexes [lsort -integer -unique $indexes]
    if {[llength $indexes] == 0} { return 0 }
    if {![info exists representation_registered($rep_key)]} {
        set representation_registered($rep_key) 1
        lappend representation_keys $rep_key
        lappend pending_representation_keys $rep_key
        set representation_kind($rep_key) $kind
        set representation_color($rep_key) $color_id
        set representation_material($rep_key) $material
        set representation_radius_tier($rep_key) $radius_tier
    } elseif {$representation_kind($rep_key) ne $kind ||
              $representation_color($rep_key) != $color_id ||
              $representation_material($rep_key) ne $material ||
              $representation_radius_tier($rep_key) ne $radius_tier} {
        error "Conflicting SQQ representation definition for '$rep_key'"
    }
    if {[info exists representation_atoms($rep_key)]} {
        set indexes [lsort -integer -unique [concat $representation_atoms($rep_key) $indexes]]
    }
    set representation_atoms($rep_key) $indexes
    return 1
}

proc ::SQQ::create_stable_representations {} {
    variable molid
    variable representation_color
    variable representation_keys
    variable representation_kind
    variable representation_material
    variable representation_name_by_key
    variable representation_radius_tier
    variable representation_style_snapshot
    variable representation_style_template
    variable pending_representation_keys
    ::SQQ::refresh_representation_style_templates
    set cage_tiers [::SQQ::active_cage_radius_tiers]
    foreach rep_key $pending_representation_keys {
        set kind $representation_kind($rep_key)
        set style_class [::SQQ::representation_style_class $rep_key $kind \
            $representation_radius_tier($rep_key)]
        if {[info exists representation_style_template($style_class)]} {
            lassign $representation_style_template($style_class) \
                representation_style color_style material_style
            mol representation {*}$representation_style
            mol color {*}$color_style
            mol material $material_style
        } else {
            if {$kind eq "cage"} {
                set radius [::SQQ::cage_layer_radius $representation_radius_tier($rep_key) $cage_tiers]
                mol representation DynamicBonds 3.5 $radius 12.0
            } elseif {$kind eq "bonds"} {
                mol representation DynamicBonds 3.5 0.125 12.0
            } elseif {$kind eq "guest"} {
                mol representation CPK 1.0 0.3 12.0 12.0
            } elseif {$kind eq "component"} {
                mol representation CPK 0.7 0.2 12.0 12.0
            } else {
                error "Unknown SQQ representation kind '$kind'"
            }
            mol color ColorID $representation_color($rep_key)
            mol material $representation_material($rep_key)
        }
        mol selection "none"
        mol addrep $molid
        set rep [expr {[molinfo $molid get numreps] - 1}]
        ::SQQ::track_representation $rep
        set representation_name_by_key($rep_key) [mol repname $molid $rep]
        set style [::SQQ::current_representation_style $rep]
        if {[llength $style] == 3} { set representation_style_snapshot($rep_key) $style }
    }
    set pending_representation_keys {}
}

proc ::SQQ::update_representation_selections {frame} {
    variable molid
    variable representation_atoms
    variable representation_keys
    variable representation_name_by_key
    foreach rep_key $representation_keys {
        if {![info exists representation_name_by_key($rep_key)]} { continue }
        set rep_name $representation_name_by_key($rep_key)
        if {[catch {mol repindex $molid $rep_name} rep] || $rep < 0} { continue }
        if {[info exists representation_atoms($rep_key)]} {
            set indexes $representation_atoms($rep_key)
        } else {
            set indexes {}
        }
        set selection [expr {[llength $indexes] == 0 ?
            "none" : "index [join $indexes { }]"}]
        mol modselect $rep $molid $selection
    }
}

proc ::SQQ::clear_graphics {} {
    variable graphics_targets
    variable graphics_ids
    variable molid
    if {$molid >= 0} {
        foreach graphics_id $graphics_ids {
            catch {graphics $molid delete $graphics_id}
        }
    }
    set graphics_ids {}
    array unset graphics_targets
    array set graphics_targets {}
}

proc ::SQQ::render_overlays {frame} {
    variable cage_centers
    variable graphics_ids
    variable graphics_targets
    variable label_visible
    variable molid
    variable pick_mode
    if {!$label_visible && $pick_mode ne "center"} { return }
    set prefix "$frame,"
    foreach key [lsort -dictionary [array names cage_centers "${frame},*"]] {
        set object_id [string range $key [string length $prefix] end]
        set center $cage_centers($key)
        if {[llength $center] != 3} { continue }
        if {$pick_mode eq "center"} {
            catch {graphics $molid color yellow}
            if {![catch {graphics $molid sphere $center radius 0.35 resolution 10} graphics_id]} {
                lappend graphics_ids $graphics_id
            }
            if {![catch {graphics $molid pickpoint $center} graphics_id]} {
                lappend graphics_ids $graphics_id
                if {[catch {graphics $molid info $graphics_id} pickpoint_info] ||
                    [llength $pickpoint_info] != 3 ||
                    ![string equal -nocase [lindex $pickpoint_info 0] pickpoint] ||
                    ![string is integer -strict [lindex $pickpoint_info 2]]} {
                    error "Unexpected VMD pickpoint metadata for graphics id $graphics_id: $pickpoint_info"
                }
                set graphics_tag [lindex $pickpoint_info 2]
                set graphics_targets($graphics_tag) $object_id
            }
        }
        if {$label_visible} {
            catch {graphics $molid color black}
            if {![catch {graphics $molid text $center $object_id size 1.0 thickness 1.0} graphics_id]} {
                lappend graphics_ids $graphics_id
            }
        }
    }
}

proc ::SQQ::guest_cages_for_atom {frame atom_index} {
    variable group_atoms
    variable group_keys
    set group_key "$frame,guest-id"
    if {![info exists group_keys($group_key)]} { return {} }
    set candidates {}
    foreach object_id $group_keys($group_key) {
        set atom_key "$frame,guest-id,$object_id"
        if {[info exists group_atoms($atom_key)] &&
            [lsearch -integer -exact $group_atoms($atom_key) $atom_index] >= 0} {
            lappend candidates $object_id
        }
    }
    return [lsort -dictionary -unique $candidates]
}

proc ::SQQ::cancel_pending_pick {} {
    variable pick_after_id
    if {$pick_after_id ne ""} {
        catch {after cancel $pick_after_id}
        set pick_after_id ""
    }
}

proc ::SQQ::remove_new_atom_labels {} {
    variable atom_label_count
    if {[catch {label list Atoms} atom_labels]} { return }
    while {[llength $atom_labels] > $atom_label_count} {
        catch {label delete Atoms [expr {[llength $atom_labels] - 1}]}
        if {[catch {label list Atoms} updated_labels]} { return }
        if {[llength $updated_labels] >= [llength $atom_labels]} { return }
        set atom_labels $updated_labels
    }
}

proc ::SQQ::schedule_atom_label_cleanup {} {
    variable label_after_ids
    set pending_ids {}
    foreach callback_id $label_after_ids {
        if {![catch {after info $callback_id}]} { lappend pending_ids $callback_id }
    }
    set label_after_ids $pending_ids
    foreach delay {0 100 250 500 1000} {
        lappend label_after_ids [after $delay [list ::SQQ::remove_new_atom_labels]]
    }
}

proc ::SQQ::apply_pick_selection {expected_mode cages guest message} {
    variable pick_after_id
    variable pick_mode
    variable selected_cages
    variable selected_guest
    set pick_after_id ""
    if {$pick_mode ne $expected_mode} { return }
    ::SQQ::remove_new_atom_labels
    set selected_cages $cages
    set selected_guest $guest
    set frame [molinfo $::SQQ::molid get frame]
    ::SQQ::render_selected $frame
    display update
    if {$message ne ""} { puts $message }
}

proc ::SQQ::queue_pick_selection {expected_mode cages guest message} {
    variable pick_after_id
    if {$pick_after_id ne ""} { catch {after cancel $pick_after_id} }
    set pick_after_id [after 50 [list \
        ::SQQ::apply_pick_selection $expected_mode $cages $guest $message]]
}

proc ::SQQ::pick_atom_event {name1 name2 operation} {
    variable atom_guest
    variable molid
    variable pick_mode
    if {$pick_mode ni {center guest} || $molid < 0} { return }
    if {![info exists ::vmd_pick_atom] ||
        ![string is integer -strict $::vmd_pick_atom] || $::vmd_pick_atom < 0} {
        return
    }
    if {![info exists ::vmd_pick_mol] || $::vmd_pick_mol != $molid} { return }
    if {$pick_mode eq "center"} {
        ::SQQ::schedule_atom_label_cleanup
        return
    }
    set frame [molinfo $molid get frame]
    set atom_key "$frame,$::vmd_pick_atom"
    if {![info exists atom_guest($atom_key)]} { return }
    set identifier $atom_guest($atom_key)
    set cages [::SQQ::guest_cages_for_atom $frame $::vmd_pick_atom]
    if {[llength $cages] == 0} {
        ::SQQ::queue_pick_selection guest {} "" \
            "SQQ pick guest: $identifier has no cage membership (frame $frame)"
        return
    }
    ::SQQ::queue_pick_selection guest $cages $identifier \
        "SQQ selected guest: $identifier; cage memberships: [join $cages { }] (frame $frame)"
}

proc ::SQQ::pick_graphics_changed {name1 name2 operation} {
    variable graphics_targets
    variable molid
    variable pick_mode
    if {$pick_mode ne "center" || $molid < 0 ||
        ![info exists ::vmd_pick_graphics]} {
        return
    }
    set values $::vmd_pick_graphics
    if {[llength $values] != 4} { return }
    lassign $values picked_molid graphics_tag button shift_state
    if {![string is integer -strict $picked_molid] ||
        $picked_molid != $molid ||
        ![string is integer -strict $graphics_tag] ||
        ![info exists graphics_targets($graphics_tag)]} {
        return
    }
    set frame [molinfo $molid get frame]
    set object_id $graphics_targets($graphics_tag)
    ::SQQ::queue_pick_selection center [list $object_id] "" \
        "SQQ selected center: $object_id (frame $frame)"
}
proc ::SQQ::render_guest_pick_context {frame} {
    variable guest_atoms
    variable guest_keys
    variable guest_types
    if {![info exists guest_keys($frame)]} { return }
    array set color_atoms {}
    foreach identifier [lsort -dictionary -unique $guest_keys($frame)] {
        set key "$frame,$identifier"
        if {![info exists guest_atoms($key)] ||
            ![info exists guest_types($key)]} {
            continue
        }
        set color_id [::SQQ::color_id guest $guest_types($key)]
        foreach atom_index $guest_atoms($key) {
            lappend color_atoms($color_id) $atom_index
        }
    }
    foreach color_id [lsort -integer [array names color_atoms]] {
        set rep_key [list pick-context guest $color_id]
        ::SQQ::register_stable_representation $frame $rep_key guest \
            [lsort -integer -unique $color_atoms($color_id)] \
            $color_id Transparent
    }
}

proc ::SQQ::create_pick_highlight_representations {} {
    variable molid
    variable pick_cage_rep_name
    variable pick_guest_rep_name
    variable pick_mode
    set pick_cage_rep_name ""
    set pick_guest_rep_name ""
    if {$pick_mode eq "off"} { return }

    mol representation DynamicBonds 3.5 0.250 12.0
    mol color ColorID [::SQQ::color_value yellow]
    mol selection "none"
    mol material Opaque
    mol addrep $molid
    set rep [expr {[molinfo $molid get numreps] - 1}]
    ::SQQ::track_representation $rep
    set pick_cage_rep_name [mol repname $molid $rep]

    if {$pick_mode eq "guest"} {
        mol representation CPK 1.0 0.3 12.0 12.0
        mol color ColorID [::SQQ::color_value orange]
        mol selection "none"
        mol material Opaque
        mol addrep $molid
        set rep [expr {[molinfo $molid get numreps] - 1}]
        ::SQQ::track_representation $rep
        set pick_guest_rep_name [mol repname $molid $rep]
    }
}

proc ::SQQ::set_pick_rep_selection {rep_name indexes} {
    variable molid
    if {$rep_name eq ""} { return }
    if {[catch {mol repindex $molid $rep_name} rep] || $rep < 0} { return }
    set selection [expr {[llength $indexes] == 0 ?
        "none" : "index [join [lsort -integer -unique $indexes] { }]"}]
    mol modselect $rep $molid $selection
}

proc ::SQQ::render_selected {frame} {
    variable group_atoms
    variable guest_atoms
    variable pick_cage_rep_name
    variable pick_guest_rep_name
    variable pick_mode
    variable selected_cages
    variable selected_guest
    if {$pick_mode eq "off"} { return }

    set cage_indexes {}
    foreach object_id $selected_cages {
        set atom_key "$frame,cage-id,$object_id"
        if {[info exists group_atoms($atom_key)]} {
            lappend cage_indexes {*}$group_atoms($atom_key)
        }
    }
    ::SQQ::set_pick_rep_selection $pick_cage_rep_name $cage_indexes

    set guest_indexes {}
    if {$pick_mode eq "guest" && $selected_guest ne ""} {
        set key "$frame,$selected_guest"
        if {[info exists guest_atoms($key)]} { set guest_indexes $guest_atoms($key) }
    }
    ::SQQ::set_pick_rep_selection $pick_guest_rep_name $guest_indexes
}

proc ::SQQ::save_target {} {
    variable cage_ids
    variable gro_path
    variable molid
    variable pick_mode
    variable selected_cages
    variable selected_guest
    if {[llength $selected_cages] == 0} {
        error "No SQQ cage is selected; use sqq pick center or sqq pick guest first"
    }
    set frame [molinfo $molid get frame]
    set targets {}
    foreach object_id $selected_cages {
        set key "$frame,$object_id"
        if {![info exists cage_ids($key)]} {
            error "Selected SQQ cage has no saved cage ID: $object_id"
        }
        lappend targets $cage_ids($key)
    }
    set targets [lsort -dictionary -unique $targets]
    set target_path [file join [file dirname $gro_path] "sqq_target.txt"]
    set temporary "${target_path}.tmp-[pid]"
    set handle [open $temporary w]
    fconfigure $handle -encoding ascii -translation lf
    puts $handle "target\t[join $targets ,]"
    puts $handle "frame\t$frame"
    puts $handle "pick_mode\t$pick_mode"
    puts $handle "guest\t$selected_guest"
    close $handle
    file rename -force $temporary $target_path
    catch {puts "SQQ target saved: [join $targets ,] -> $target_path"}
    return $target_path
}
proc ::SQQ::announce_graph_mode {frame} {
    variable graph_mode
    variable displayed_graph_mode
    set value [expr {[info exists graph_mode($frame)] ? $graph_mode($frame) : "unknown"}]
    if {$value ne $displayed_graph_mode} {
        puts "SQQ graph: $value"
        set displayed_graph_mode $value
    }
}

proc ::SQQ::ordered_active_families {} {
    variable active_families
    set ordered {}
    foreach family {component phase cluster domain cage guest} {
        if {$family in $active_families} { lappend ordered $family }
    }
    return $ordered
}

proc ::SQQ::render_family {frame family targets} {
    variable group_atoms
    set representation_count 0
    set material [::SQQ::base_material]
    if {$family in {cage guest}} {
        array set explicit_ids {}
        foreach target $targets {
            lassign $target source key
            if {$source eq "${family}-id" || $source eq "${family}-track"} {
                set explicit_ids($key) 1
            }
        }
        array set layer_atoms {}
        set layer_keys {}
        foreach item [::SQQ::expanded_targets $frame $family $targets] {
            lassign $item source key
            set atom_key [::SQQ::group_atom_key $frame $source $key]
            if {![info exists group_atoms($atom_key)]} { continue }
            lassign [::SQQ::effective_color $frame $source $key] color_id color_priority
            set explicit [info exists explicit_ids($key)]
            set layer_key [::SQQ::object_render_key $frame $key $color_priority $color_id $explicit]
            if {![info exists layer_atoms($layer_key)]} {
                lappend layer_keys $layer_key
                set layer_atoms($layer_key) {}
            }
            foreach atom_index $group_atoms($atom_key) { lappend layer_atoms($layer_key) $atom_index }
        }
        set layer_keys [lsort -command ::SQQ::compare_object_render_keys $layer_keys]
        if {$family eq "cage"} {
            foreach layer_key $layer_keys {
                set indexes [lsort -integer -unique $layer_atoms($layer_key)]
                set rep_key [list main $family $layer_key]
                incr representation_count [::SQQ::register_stable_representation \
                    $frame $rep_key cage $indexes [lindex $layer_key 6] $material \
                    [::SQQ::cage_radius_tier $layer_key]]
            }
        } else {
            foreach layer_key $layer_keys {
                set indexes [lsort -integer -unique $layer_atoms($layer_key)]
                set rep_key [list main $family $layer_key]
                incr representation_count [::SQQ::register_stable_representation \
                    $frame $rep_key guest $indexes [lindex $layer_key 6] $material]
            }
        }
    } elseif {$family eq "component"} {
        array set color_atoms {}
        set render_keys {}
        foreach item [::SQQ::expanded_targets $frame $family $targets] {
            lassign $item source key
            set atom_key [::SQQ::group_atom_key $frame $source $key]
            if {![info exists group_atoms($atom_key)]} { continue }
            lassign [::SQQ::effective_color $frame $source $key] color_id priority
            set render_key "$priority,$color_id"
            if {![info exists color_atoms($render_key)]} {
                lappend render_keys [list $priority $color_id]
            }
            foreach atom_index $group_atoms($atom_key) {
                lappend color_atoms($render_key) $atom_index
            }
        }
        foreach render_key [lsort -command ::SQQ::compare_render_keys $render_keys] {
            lassign $render_key priority color_id
            set indexes [lsort -integer -unique $color_atoms($priority,$color_id)]
            set rep_key [list main $family $render_key]
            incr representation_count [::SQQ::register_stable_representation \
                $frame $rep_key component $indexes $color_id $material]
        }
    } else {
        array set color_atoms {}
        set render_keys {}
        foreach item [::SQQ::expanded_targets $frame $family $targets] {
            lassign $item source key
            set atom_key [::SQQ::group_atom_key $frame $source $key]
            if {![info exists group_atoms($atom_key)]} { continue }
            lassign [::SQQ::effective_color $frame $source $key] color_id priority
            set render_key "$priority,$color_id"
            if {![info exists color_atoms($render_key)]} { lappend render_keys [list $priority $color_id] }
            foreach atom_index $group_atoms($atom_key) { lappend color_atoms($render_key) $atom_index }
        }
        foreach render_key [lsort -command ::SQQ::compare_render_keys $render_keys] {
            lassign $render_key priority color_id
            set indexes [lsort -integer -unique $color_atoms($priority,$color_id)]
            set rep_key [list main $family $render_key]
            incr representation_count [::SQQ::register_stable_representation \
                $frame $rep_key bonds $indexes $color_id $material]
        }
    }
    return $representation_count
}

proc ::SQQ::render_current {{announce 0}} {
    ::SQQ::cancel_pending_render
    variable molid
    variable active_families
    variable active_targets
    variable pick_mode
    if {$molid < 0 || $molid ni [molinfo list]} { return }
    set frame [molinfo $molid get frame]
    ::SQQ::announce_graph_mode $frame
    ::SQQ::clear_representations
    ::SQQ::clear_graphics
    ::SQQ::initialize_component_frame $frame
    variable representation_atoms
    array unset representation_atoms
    array set representation_atoms {}
    array set current_counts {}
    foreach family [::SQQ::ordered_active_families] {
        if {$pick_mode eq "guest" && $family eq "guest"} { continue }
        set targets $active_targets($family)
        set current_counts($family) [::SQQ::render_family $frame $family $targets]
    }
    if {$pick_mode eq "guest"} {
        ::SQQ::render_guest_pick_context $frame
    }
    ::SQQ::create_stable_representations
    if {$announce} {
        foreach family [::SQQ::ordered_active_families] {
            if {$pick_mode eq "guest" && $family eq "guest"} { continue }
            set targets $active_targets($family)
            set labels {}
            foreach target $targets { lappend labels [::SQQ::target_label $family $target] }
            set representation_count [expr {[info exists current_counts($family)] ?
                $current_counts($family) : 0}]
            if {$representation_count == 0} {
                puts "SQQ show $family: no memberships for [join $labels { }] in frame $frame"
            } else {
                puts "SQQ show $family: [join $labels { }] (frame $frame)"
            }
        }
    }
    ::SQQ::create_pick_highlight_representations
    ::SQQ::update_representation_selections $frame
    ::SQQ::render_selected $frame
    ::SQQ::render_overlays $frame
    display update
}

proc ::SQQ::update_current {} {
    ::SQQ::cancel_pending_render
    variable molid
    if {$molid < 0 || $molid ni [molinfo list]} { return }
    set frame [molinfo $molid get frame]
    ::SQQ::announce_graph_mode $frame
    ::SQQ::initialize_component_frame $frame
    variable active_targets
    variable pick_mode
    variable representation_atoms
    array unset representation_atoms
    array set representation_atoms {}
    foreach family [::SQQ::ordered_active_families] {
        if {$pick_mode eq "guest" && $family eq "guest"} { continue }
        ::SQQ::render_family $frame $family $active_targets($family)
    }
    if {$pick_mode eq "guest"} {
        ::SQQ::render_guest_pick_context $frame
    }
    ::SQQ::create_stable_representations
    ::SQQ::update_representation_selections $frame
    ::SQQ::clear_graphics
    ::SQQ::render_selected $frame
    ::SQQ::render_overlays $frame
    display update
}

proc ::SQQ::cancel_pending_render {} {
    variable frame_after_id
    if {$frame_after_id ne ""} {
        catch {after cancel $frame_after_id}
        set frame_after_id ""
    }
}

proc ::SQQ::render_pending {} {
    variable frame_after_id
    set frame_after_id ""
    ::SQQ::update_current
}

proc ::SQQ::frame_changed {name1 name2 operation} {
    ::SQQ::cancel_pending_pick
    variable molid
    variable frame_after_id
    variable selected_cages
    variable selected_guest
    if {$name2 ne "$molid"} { return }
    set selected_cages {}
    set selected_guest ""
    if {$frame_after_id ne ""} { catch {after cancel $frame_after_id} }
    set frame_after_id [after idle [list ::SQQ::render_pending]]
}

proc ::SQQ::startup_banner {} {
__SQQ_BANNER_BODY__
    puts ""
}

proc ::SQQ::startup_help {} {
    puts "Default view : cage all (opaque)"
    puts "Show mode    : additive"
    puts ""
    puts "Commands:"
    puts "  sqq show <family> <target...> ?<family> <target...> ...?"
    puts "  sqq show label ?on|off?"
    puts "  sqq color <family> <target...> <color>"
    puts "  sqq pick center|guest|off"
    puts "  sqq target save"
    puts "  sqq clear"
    puts "  sqq -h"
}

proc ::SQQ::help {} {
__SQQ_HELP_BODY__
}

proc sqq {{command help} args} {
    set normalized_command [string tolower $command]
    if {$normalized_command ni {help -h --help}} {
        ::SQQ::require_renderer_ready
    }
    switch -- $normalized_command {
        show {
            if {[llength $args] >= 1 &&
                [string tolower [lindex $args 0]] in {label lable}} {
                ::SQQ::set_label [lrange $args 1 end]
            } else {
                if {[llength $args] < 2} {
                    error "Usage: sqq show <family> <target...> ?<family> <target...> ...?"
                }
                ::SQQ::set_show $args
            }
        }
        color {
            if {[llength $args] < 3} {
                error "Usage: sqq color <family> <target> ?target ...? <VMD-color|ColorID|default>"
            }
            ::SQQ::set_colors [lindex $args 0] [lrange $args 1 end-1] [lindex $args end]
        }
        pick {
            if {[llength $args] != 1} { error "Usage: sqq pick center|guest|off" }
            ::SQQ::set_pick_mode [lindex $args 0]
        }
        target {
            if {[llength $args] != 1 || ![string equal -nocase [lindex $args 0] save]} {
                error "Usage: sqq target save"
            }
            ::SQQ::save_target
        }
        clear {
            if {[llength $args] != 0} { error "Usage: sqq clear" }
            ::SQQ::reset_show 1
        }
        help - -h - --help {
            if {[llength $args] != 0} { error "Usage: sqq help" }
            ::SQQ::help
        }
        default { error "Unknown SQQ command '$command'; use show, color, pick, target, clear, or help" }
    }
}

set script_dir [file dirname [file normalize [info script]]]
set ::SQQ::gro_path [file join $script_dir __SQQ_GRO_FILENAME__]
set ::SQQ::xtc_path [file join $script_dir __SQQ_XTC_FILENAME__]
set ::SQQ::membership_path [file join $script_dir __SQQ_MEMBERSHIP_FILENAME__]
set missing_render_file ""
foreach {label path} [list GRO $::SQQ::gro_path XTC $::SQQ::xtc_path membership $::SQQ::membership_path] {
    if {![file isfile $path]} {
        set missing_render_file "SQQ $label file not found: $path"
        break
    }
}
::SQQ::startup_banner
catch {color Display Background white}
if {$missing_render_file ne ""} {
    ::SQQ::renderer_fail $missing_render_file
    error "SQQ renderer failed: $missing_render_file"
} else {
    puts "SQQ VMD Renderer: loading files, please wait..."
    if {[catch {::SQQ::load_renderer} message options]} {
        ::SQQ::renderer_fail $message
        return -options $options "SQQ renderer failed: $message"
    }
}
"""

__all__ = ["SQQ_CAGE_TCL"]
